import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import scipy.signal

from sklearn import linear_model

class Subject:

    def __init__(self, importpath, importpath_evt, importpath_csv):
        datafile = sp.io.loadmat(importpath)
        self.name   = str(np.squeeze(datafile['name']))
        self.srate  = int(np.squeeze(datafile['srate']))
        self.data   = np.squeeze(datafile['data'])
        self.nbchan = len(self.data)
        self.events = {}
        self.events['df'] = pd.read_csv(importpath_evt, sep='\t')
        self._construct_event_hierarchy()


    def _construct_event_hierarchy(self):
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


    def modify_trial_length(self, new_trial_length):
        """
        Reduces or increases trial lengths and modifies clean segment
        markers to reflect new trial lengths.
        """
        df = self.events['df'].copy()
        for i in range(df.shape[0]):
            if df.iloc[i].Type in ['01', '11']:
                for j in range(i+1, df.shape[0]):
                    if df.iloc[j].Type in ['02', '12']:
                        break
                df.at[j, 'Latency'] = df.iloc[i].Latency + new_trial_length
        self.events['df'] = df.copy()
        self._update_event_hierarchy()


    def get_windows(self, chan, seg_type, nperwindow=512*2, noverlap=512):
        """ Grabs windows of data of size nperwindow with overlap noverlap.
        TODO: Update. This is technically not a correct implementation of Welch's method,
        but it works if noverlap is 50 percent of the window length.
        Arguments
            chan:       Integer, channel number from which to extract windows.
            seg_type:   String, specifies what segments to use. Possible options are
                        'C' or 'O'.
            nperwindow: Time-points to use per window. Default value, provided sampling rate
                        is 512 Hz, is 2 seconds.
            noverlap:   Overlap of windows. Default is 50%.
        """
        wins = []
        for event in self.events[seg_type]:
            event_length = event[1] - event[0]
            if event_length >= nperwindow:
                current_idx = 0
                while current_idx + nperwindow <= event_length:
                    wins.append(self.data[chan][current_idx : current_idx+nperwindow])
                    current_idx += noverlap
        return wins


    def welch(self, windows, srate):
        """
        Takes a list of data segments (each size 1xN), computes each segment's PSD,
        and averages them to get a final PSD.
        """
        # psds = [sp.signal.periodogram(window, srate, window='hamming')[1] for window in windows]
        psds = [sp.signal.welch(window, srate, nperseg=len(window), window='hamming')[1] for window in windows]
        return np.mean(psds, axis=0)


    def remove_freq_buffer(self, data, lofreq, hifreq):
        """
        Removes a frequency buffer from a PSD or frequency vector.
        """
        data = np.delete(data, range(lofreq*2, hifreq*2))
        return data.reshape(len(data), 1)


    def compute_ch_psds(self, match_OA_protocol=True, manual_win_extraction=False, 
                                            nwins_upperlimit=-1, rand_wins=False):
        """ Returns subj data structure with calculated PSDs and subject information.
        Arguments:
            import_path:          String, path to .mat files
            import_path_csv:      String, path to .csv containing subject class, sex, and
                                  age information.
            cut_recording_length: Boolean, specifies whether to cut recordings down to 7
                                  minutes.
            nwins_upperlimit:     Scalar, specifies the upperlimit of windows to extract
                                  from a given subject. A value of -1 means no upper limit.
        """

        self.psds = {}
        self.f = np.linspace(0, 256, 513)
        self.f = self.f.reshape(len(self.f), 1)
        self.f_rm_alpha = self.remove_freq_buffer(self.f, 7, 14)
        
        print('Computing PSDs: {}...'.format(self.name), end='')

        # Cut trials down to match older adults; 30 seconds
        if match_OA_protocol and self.name[0:3] == '112':
            self.modify_trial_length(512 * 30)

        for ch in range(self.nbchan):
            self.psds[ch] = {}

            if manual_win_extraction:
                eyesc_windows = self.get_windows(ch, 'eyesc')
                eyeso_windows = self.get_windows(ch, 'eyeso')

                # Remove windows if there's an upper-limit and windows
                # are selected randomly.
                if nwins_upperlimit != -1 and rand_wins:
                    while len(eyesc_windows) > nwins_upperlimit:
                        if rand_wins:
                            random.shuffle(eyesc_wins)
                            eyesc_wins.pop()
                        else:
                            eyesc_windows.pop()
                    while len(eyeso_windows) > nwins_upperlimit:
                        if rand_wins:
                            random.shuffle(eyeso_windows)
                            eyeso_windows.pop(random.randrange(eyeso_windows))
                        else:
                            eyeso_windows.pop()

                self.psds[ch]['eyesc'] = self.welch(eyesc_windows, 512)
                self.psds[ch]['eyeso'] = self.welch(eyeso_windows, 512)
            else:
                eyesc_segs = [self.data[ch][event[0] : event[1]] for event in self.events['eyesc']]
                eyeso_segs = [self.data[ch][event[0] : event[1]] for event in self.events['eyeso']]
                eyesc_psds = []
                eyeso_psds = []
                for seg in eyesc_segs:
                    if len(seg) >= 1024:
                        eyesc_psds.append(sp.signal.welch(seg, 512, nperseg=1024, noverlap=512, window='hamming')[1])
                for seg in eyeso_segs:
                    if len(seg) >= 1024:
                        eyeso_psds.append(sp.signal.welch(seg, 512, nperseg=1024, noverlap=512, window='hamming')[1])
                self.psds[ch]['eyesc'] = np.mean(eyesc_psds, axis=0)
                self.psds[ch]['eyeso'] = np.mean(eyeso_psds, axis=0)

        if manual_win_extraction:
            self.nwins_eyesc = len(eyesc_windows)
            self.nwins_eyeso = len(eyeso_windows)
        else:
            self.nwins_eyesc = 0
            self.nwins_eyeso = 0
            for event in self.events['eyesc']:
                event_length = event[1] - event[0]
                if event_length >= 1024:
                    self.nwins_eyesc += (event_length//512) - 1
            for event in self.events['eyeso']:
                event_length = event[1] - event[0]
                if event_length >= 1024:
                    self.nwins_eyeso += (event_length//512) - 1
        self.data = [] # Clear it from memory since it's no longer needed.
        print('Done.')


    def linreg_slope(self, ch, lofreq, hifreq):
        model = linear_model.LinearRegression()
        model.fit(self.f_rm_alpha[lofreq*2:hifreq*2], np.log10(self.psds[ch]['eyesc_rm_alpha'][lofreq*2:hifreq*2]))
        fit_line = model.predict(self.f_rm_alpha)
        return model.coef_[0] * (10**2), fit_line


    def ransac_slope(self, ch, lofreq, hifreq):
        """
        Robustly fits line to the PSD, using the RANSAC algorithm.
        Returns slope and fit line.
        """
        model = linear_model.RANSACRegressor(linear_model.LinearRegression())
        model.fit(self.f_rm_alpha[lofreq*2:hifreq*2], np.log10(self.psds[ch]['eyesc_rm_alpha'][lofreq*2:hifreq*2]))
        fit_line = model.predict(self.f_rm_alpha)
        return model.estimator_.coef_[0] * (10**2), fit_line


    def fit_slopes(self, regr_func_str, lofreq, hifreq):
        print('Fitting: {}...'.format(self.name), end='')
        if regr_func_str == 'ransac':
            regr_func = self.ransac_slope
        elif regr_func_str == 'linreg':
            regr_func = self.linreg_slope
        for ch in range(self.nbchan):
            self.psds[ch]['eyesc_rm_alpha'] = self.remove_freq_buffer(self.psds[ch]['eyesc'], 7, 14)
            self.psds[ch]['eyeso_rm_alpha'] = self.remove_freq_buffer(self.psds[ch]['eyeso'], 7, 14)
            self.psds[ch]['eyesc_slope'], self.psds[ch]['eyesc_fitline'] = regr_func(ch, lofreq, hifreq)
            self.psds[ch]['eyeso_slope'], self.psds[ch]['eyeso_fitline'] = regr_func(ch, lofreq, hifreq)
        print('Done.')

