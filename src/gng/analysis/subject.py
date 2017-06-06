import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import scipy.signal

from sklearn import linear_model
from events import rm_intertrial_segs

class Subject:

    def __init__(self, importpath, importpath_evt, group='', age=0, sex=0):
        datafile = sp.io.loadmat(importpath)
        self.name   = str(np.squeeze(datafile['name']))
        self.group  = group
        self.age    = age
        self.sex    = sex
        self.srate  = int(np.squeeze(datafile['srate']))
        self.data   = np.squeeze(datafile['data'])
        self.nbchan = len(self.data)
        self.events = {}
        self.events['df'] = pd.read_csv(importpath_evt, sep='\t')
        self._construct_event_hierarchy()


    def _construct_event_hierarchy(self, clean_space='all'):
        """
        Constructs trial hierarchy.
        1. Based on the value of clean_space, construct trial hierarchy.
        2. If 'all':
            take entire intertrial space
        """






        self.events['df'].Type = self.events['df'].Type.map(
            lambda x: x[0:2][::-1] if len(x) == 3 else x
        )
        self.events['types'] = set(map(lambda x: x[0], self.events['df'].Type))
        for event_type in self.events['types']:
            if event_type == 'C':
                self._organize_events_into_df('eyesc', event_type)
            elif event_type == 'O':
                self._organize_events_into_df('eyeso', event_type)
            elif event_type == '0':
                self._organize_events_into_df('trials_eyesc', event_type)
            elif event_type == '1':
                self._organize_events_into_df('trials_eyeso', event_type)


    def _organize_events_into_df(self, event_label, event_type):
        self.events[event_label] = []
        idx = list(map(lambda x: True if x[0] == event_type\
                                      else False, self.events['df'].Type))
        events = self.events['df'][idx]
        for i in range(events.shape[0]):
            if events.iloc[i].Type == event_type + '1':
                self.events[event_label].append(
                    [events.iloc[i].Latency, events.iloc[i+1].Latency]
                )


    def _update_event_hierarchy(self):
        """
        Checks to make sure that clean segments fall inside of trials.
        For example, run this after running modify_trial_length().
        """
        self.events['df'] = rm_intertrial_segs(self.events['df'])
        self._construct_event_hierarchy()


    def modify_trial_length(self, lower_bound, upper_bound):
        """
        Reduces or increases trial lengths and modifies clean segment
        markers to reflect new trial lengths.
        Arguments
            lower_bound : float (seconds)
                Specifies lower bound for new trial length
            upper_bound : float (seconds)
                Specifies upper bound for new trial length
        For example, limiting 60-second trials down to 30-second ones
        would be done like so:
            >> s = Subject('1121181181.mat', '1121181181.evt')
            >> s.modify_trial_length(0, 30)
        """
        df = self.events['df'].copy()
        for i in range(df.shape[0]):
            if df.iloc[i].Type in ['01', '11']:
                for j in range(i+1, df.shape[0]):
                    if df.iloc[j].Type in ['02', '12']:
                        break
                df.at[j, 'Latency'] = df.iloc[i].Latency + upper_bound*512
                df.at[i, 'Latency'] = df.iloc[i].Latency + lower_bound*512
        self.events['df'] = df.copy()
        self._update_event_hierarchy()


    def get_windows(self, chan, nperwindow=512*2, noverlap=512):
        """ Grabs windows of data of size nperwindow with overlap noverlap.
        TODO: Update. This is technically not a correct implementation of Welch's method,
        but it works if noverlap is 50 percent of the window length.
        Arguments
            chan:       Integer, channel number from which to extract windows.
            nperwindow: Time-points to use per window. Default value, provided sampling rate
                        is 512 Hz, is 2 seconds.
            noverlap:   Overlap of windows. Default is 50%.
        """
        segs = []
        for event in self.events[seg_type]:
            segs.append(self.data[chan][event[0] : event[1]])

        wins = []
        for seg in segs:
            current_idx = 0
            while current_idx + nperwindow <= len(seg):
                wins.append(seg[current_idx : current_idx+nperwindow])
                current_idx += noverlap
        return wins


    def welch(self, windows, srate):
        """
        Takes a list of data segments (each size 1xN), computes each segment's PSD,
        and averages them to get a final PSD.
        """
        psds = [sp.signal.periodogram(w, srate, window='hamming')[1] for w in windows]
        return np.mean(psds, axis=0)


    def remove_freq_buffer(self, data, lofreq, hifreq):
        """
        Removes a frequency buffer from a PSD or frequency vector.
        """
        data = np.delete(data, range(lofreq*2, hifreq*2))
        return data.reshape(len(data), 1)


    def compute_ch_psds(self, nwins_upperlimit=0):
        """
        Returns subj data structure with calculated PSDS and subject
        information.
        Arguments:
            nwins_upperlimit : int
                Upper limit on number of windows we use to compute the
                PSD. Default is 0, which means no upper limit.
        """
        self.psds = {}
        self.f = np.linspace(0, 256, 513)
        self.f = self.f.reshape(len(self.f), 1)
        self.f_rm_alpha = self.remove_freq_buffer(self.f, 7, 14)

        for ch in range(self.nbchan):
            self.psds[ch] = {}
            windows = self.get_windows(ch)
            if nwins_upperlimit:
                while len(eyesc_windows) > nwins_upperlimit:
                    random.shuffle(eyesc_windows)
                    random.shuffle(eyeso_windows)
                    eyesc_windows.pop()
                    eyeso_windows.pop()
            self.psds[ch]['eyesc'] = self.welch(eyesc_windows, self.srate)
            self.psds[ch]['eyeso'] = self.welch(eyeso_windows, self.srate)

        self.nwins_eyesc = len(eyesc_windows)
        self.nwins_eyeso = len(eyeso_windows)
        self.data = [] # Clear it from memory since it's no longer needed.


    def linreg_slope(self, f, psd, lofreq, hifreq):
        """
        Fits line to the PSD, using simple linear regression.
        Returns slope and fit line.
        """
        model = linear_model.LinearRegression()
        model.fit(f[lofreq*2:hifreq*2], np.log10(psd[lofreq*2:hifreq*2]))
        fit_line = model.predict(f)
        return model.coef_[0] * (10**2), fit_line


    def ransac_slope(self, f, psd, lofreq, hifreq):
        """
        Robustly fits line to the PSD, using the RANSAC algorithm.
        Returns slope and fit line.
        """
        model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
        model_ransac.fit(f[lofreq*2:hifreq*2], np.log10(psd[lofreq*2:hifreq*2]))
        fit_line = model_ransac.predict(f)
        return model_ransac.estimator_.coef_[0] * (10**2), fit_line


    def fit_slopes(self, regr_func_str='ransac',
                    buffer_lofreq=7, buffer_hifreq=14,
                    fitting_lofreq=2, fitting_hifreq=24):
        if regr_func_str == 'ransac':
            regr_func = self.ransac_slope
        elif regr_func_str == 'linreg':
            regr_func = self.linreg_slope
        for ch in range(self.nbchan):
            self.psds[ch]['eyesc_rm_alpha'] = self.remove_freq_buffer(self.psds[ch]['eyesc'], buffer_lofreq, buffer_hifreq)
            self.psds[ch]['eyeso_rm_alpha'] = self.remove_freq_buffer(self.psds[ch]['eyeso'], buffer_lofreq, buffer_hifreq)
            eyesc_slope, eyesc_fitline = regr_func(self.f_rm_alpha, self.psds[ch]['eyesc_rm_alpha'], fitting_lofreq, fitting_hifreq)
            eyeso_slope, eyeso_fitline = regr_func(self.f_rm_alpha, self.psds[ch]['eyeso_rm_alpha'], fitting_lofreq, fitting_hifreq)
            self.psds[ch]['eyesc_slope'], self.psds[ch]['eyesc_fitline'] = eyesc_slope, eyesc_fitline
            self.psds[ch]['eyeso_slope'], self.psds[ch]['eyeso_fitline'] = eyeso_slope, eyeso_fitline
