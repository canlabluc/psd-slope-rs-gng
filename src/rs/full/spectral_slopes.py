"""
Computes PSDs and fits line to specified frequency range in PSD for
BESA source models.

Change parameters below in order to run different datasets.
"""

import os
import glob
import datetime
import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import scipy.signal
from sklearn import linear_model

#### Parameters ####
# `recompute_psds`: `True` or `False`, for recomputing subject PSDs or loading previous results.
# `psd_buffer_lofreq`: Scalar, specifies the lower bound for the PSD buffer that we exclude.
# `psd_buffer_hifreq`: Scalar, specifies the upper bound for the PSD buffer that we exclude.
# `fitting_func`: `'linreg'` or `'ransac'`, specifies which function to use for fitting. `'linreg'` is simple linear regression. `'ransac'` is RANSAC, a robust fitting method that ignores outliers.
# `fitting_lofreq`: Scalar, specifies the lower bound for the PSD fitting range.
# `fitting_hifreq`: Scalar, specifies the upper bound for the PSD fitting range.
# `import_dir`: String specifying the directory to import results to.
# `export_dir`: String specifying the directory to export results to.
################################################################################

montage = 'ventral'
recompute_psds = True
psd_buffer_lofreq = 7
psd_buffer_hifreq = 14
fitting_func = 'ransac'
fitting_lofreq = 2
fitting_hifreq = 24
nwins_upperlimit = 100
import_dir = '/Users/jorge/Drive/research/_psd-slope/data/rs/full/source-ventral/MagEvtFiltCAR-mat/'
# export_dir = '/Users/jorge/Drive/research/_psd-slope/data/runs/'
export_dir = '/Users/jorge/'

################################################################################

#### Functions ####

def get_filelist(import_path):
    matfiles = []
    for root, dirs, files in os.walk(import_path):
        matfiles += glob.glob(os.path.join(root, '*.mat'))
    return matfiles

def import_subject(subj, i, import_path):
    """ Imports a single subject and adds them to the subj data structure.
    Inputs:
        - subj: Dictionary, contains all subject information. Constructed by
        compute_subject_psds()
        - i: Integer, specifies the ID of the current subject being imported.
        import_path: String, absolute path to directory containing .mat files.
    """
    subj[i] = {}
    datafile = sp.io.loadmat(import_path)
    subj[i]['name'] = str(np.squeeze(datafile['name']))
    subj[i]['srate'] = int(np.squeeze(datafile['srate']))
    subj[i]['events'] = []
    for event in np.squeeze(datafile['evts']):
        subj[i]['events'].append([event[0][0], event[1][0][0], event[2][0][0]])
    subj[i]['data'] = np.squeeze(datafile['data'])
    subj[i]['nbchan'] = len(subj[i]['data'])
    return subj

def get_num_extractable_windows(events, print_info):
    """
    Helper function that returns the total number of windows and seconds
    able to be extracted from a specified resting state recording.
    Inputs:
        - events: List of event codes, each event in the format: [port_code, time_point],
        such as ['C1', 10249]
        - port_code: String specifying the events to extract. For example, 'C1' or 'O1'
        for eyes-closed and eyes-open resting state data, respectively.
        - print_info: Boolean specifying whether to print info to the terminal.

    get_num_extractable_windows assumes that we're looking for windows of 2-second lengths,
    with 50% overlap. Thus, 3 seconds of extractable data provides us with 2 windows.
    """
    for event in events:
        # If we can extract at least two seconds...
        if (event[1] - event[0]) >= 1024:
            points = event[1] - event[0]
            secs   = (event[1] - event[0])//512
            nwin   = (event[1] - event[0])//512 - 1 # Assuming 50% window overlap.
            total_wins += nwin
            total_secs += secs
            if print_info == True:
                print('Event {}:\t{} points, {} seconds, {} windows'.format(event, points, secs, nwin))
    if print_info == True:
        print('Total windows able to be extracted: ', total_wins)
        print('Total seconds able to be extracted: ', total_secs)
    return total_wins, total_secs

def get_windows(data, events, nperwindow=512*2, noverlap=512):
    """ Grabs windows of data of type port_code using events information.
    Arguments:
        data:       A channel of data from which to extract windows.
        events:     List of time segments from which to extract data, with a time
                    segment in the format: [start_timepoint, ending_timepoint]
        port_code:  String specifying the events to extract. For example, 'C1' or 'O1'
                    for eyes-closed and eyes-open resting state data, respectively.
        nperwindow: Time-points to use per window. Default value, provided sampling rate
                    is 512 Hz, is 2 seconds.
        noverlap:   Overlap of windows. Default is 50%.
    """
    windows = []
    for event in events:
        if event[1]-event[0] >= nperwindow:
            nwindows = (event[1] - event[0])//noverlap - 1
            for i in range(nwindows):
                windows.append(data[event[0] + noverlap*i : event[0] + noverlap*i + nperwindow])
    return windows

def welch(windows, srate):
    """
    Takes a list of data segments (each size 1xN), computes each segment's PSD,
    and averages them to get a final PSD.
    """
    psds = [sp.signal.welch(window, srate, nperseg=len(window), window='hamming')[1] for window in windows]
    return np.mean(psds, axis=0)

def remove_freq_buffer(data, lofreq, hifreq):
    """
    Removes a frequency buffer from a PSD or frequency vector.
    """
    data = np.delete(data, range(lofreq*2, hifreq*2))
    return data.reshape(len(data), 1)

def compute_subject_psds(import_path, import_path_csv):
    """ Returns subj data structure with calculated PSDs and subject information.
    Arguments:
        import_path:     String, path to .mat files
        import_path_csv: String, path to .csv containing subject class, sex, and
                         age information.
    """
    matfiles = get_filelist(import_path)
    df = pd.read_csv(import_path_csv)
    df.SUBJECT = df.SUBJECT.astype(str)

    # Making sure that all subjects are available in the .csv
    subjects = set(map(lambda x: x.split('/')[-1][:-4], matfiles))
    missing = subjects - set(df.SUBJECT)
    if len(missing) != 0:
        for s in missing:
            print('ERROR: Specified csv does not contain information for subject {}.\n'.format(s))
        raise Exception('Missing subject information from csv. Either remove subject file from\n'+
                        'processing pipeline or add subject information to csv file.')

    subj = {}
    subj['nbsubj'] = len(matfiles)
    subj['f'] = np.linspace(0, 256, 513)
    subj['f'] = subj['f'].reshape(len(subj['f']), 1)
    subj['f_rm_alpha'] = remove_freq_buffer(subj['f'], 7, 14)
    for i in range(len(matfiles)):
        print('Processing: {}... '.format(matfiles[i].split('/')[-1]), end='')

        subj = import_subject(subj, i, matfiles[i])
        subj[i]['age']   = df[df.SUBJECT == subj[i]['name']].AGE.values[0]
        subj[i]['class'] = df[df.SUBJECT == subj[i]['name']].CLASS.values[0]
        subj[i]['sex']   = df[df.SUBJECT == subj[i]['name']].SEX.values[0]

        # Reorganize events into two separate lists: Eyes closed and eyes open.
        subj[i]['events_eyesc'] = [[subj[i]['events'][j][1], subj[i]['events'][j+1][1]] for j in range(len(subj[i]['events'])) if subj[i]['events'][j][0] == 'C1']
        subj[i]['events_eyeso'] = [[subj[i]['events'][j][1], subj[i]['events'][j+1][1]] for j in range(len(subj[i]['events'])) if subj[i]['events'][j][0] == 'O1']

        for ch in range(subj[i]['nbchan']):
            subj[i][ch] = {}
            eyesC_windows = get_windows(subj[i]['data'][ch], subj[i]['events_eyesc'])
            eyesO_windows = get_windows(subj[i]['data'][ch], subj[i]['events_eyeso'])
            # Discard windows from the back of the recording if the subject has more than 100.
            while len(eyesC_windows) > nwins_upperlimit:
                eyesC_windows.pop()
            while len(eyesO_windows) > nwins_upperlimit:
                eyesO_windows.pop()
            subj[i][ch]['eyesC_psd'] = welch(eyesC_windows, 512)
            subj[i][ch]['eyesO_psd'] = welch(eyesO_windows, 512)
            subj[i][ch]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['eyesC_psd'], 7, 14)
            subj[i][ch]['eyesO_psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['eyesO_psd'], 7, 14)
        subj[i]['eyesc_nwindows'] = len(eyesC_windows)
        subj[i]['eyeso_nwindows'] = len(eyesO_windows)
        subj[i]['data'] = np.nan # No longer needed, so clear it from memory
        subj[i]['eyesC_psd'] = np.mean([subj[i][ch]['eyesC_psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['eyesO_psd'] = np.mean([subj[i][ch]['eyesO_psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i]['eyesC_psd'], 7, 14)
        subj[i]['eyesO_psd_rm_alpha'] = remove_freq_buffer(subj[i]['eyesO_psd'], 7, 14)
        print('Done.')
    return subj

def linreg_slope(f, psd, lofreq, hifreq):
    """
    Fits line to the PSD, using simple linear regression.
    Returns slope and fit line.
    """
    model = linear_model.LinearRegression()
    model.fit(f[lofreq*2:hifreq*2], np.log10(psd[lofreq*2:hifreq*2]))
    fit_line = model.predict(f)
    return model.coef_[0] * (10**2), fit_line

def ransac_slope(f, psd, lofreq, hifreq):
    """
    Robustly fits line to the PSD, using the RANSAC algorithm.
    Returns slope and fit line.
    """
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
    model_ransac.fit(f[lofreq*2:hifreq*2], np.log10(psd[lofreq*2:hifreq*2]))
    fit_line = model_ransac.predict(f)
    return model_ransac.estimator_.coef_[0] * (10**2), fit_line

def fit_slopes(subj, regr_func, lofreq, hifreq):
    """
    Takes subj data structure and fits slopes to each subject's PSDs and mean
    PSD, using regr_func and fitting to datapoints between lofreq and hifreq.
    """
    for i in range(subj['nbsubj']):
        # Per-subject PSD average fitting
        subj[i]['eyesC_slope'], subj[i]['eyesC_fitline'] = regr_func(subj['f'], subj[i]['eyesC_psd'], lofreq, hifreq)
        subj[i]['eyesO_slope'], subj[i]['eyesO_fitline'] = regr_func(subj['f'], subj[i]['eyesO_psd'], lofreq, hifreq)
        for ch in range(subj[i]['nbchan']):
            # Per-channel PSD fitting
            subj[i][ch]['eyesC_slope'], subj[i][ch]['eyesC_fitline'] = regr_func(subj['f'], subj[i][ch]['eyesC_psd_rm_alpha'], lofreq, hifreq)
            subj[i][ch]['eyesO_slope'], subj[i][ch]['eyesO_fitline'] = regr_func(subj['f'], subj[i][ch]['eyesO_psd_rm_alpha'], lofreq, hifreq)
    return subj

def get_subject_slopes(subj, ch, slope_type):
    """ Returns list of slopes for specified channel of slope_type.
    Arguments:
        subj: The subj data structure.
        ch:   Scalar, channel for which to get list of subject slopes.
        slope_type: String, e.g., 'eyesO_slope' or 'eyesC_slope'
    """
    if ch == -1: # Slope of PSD grand average
        return [subj[i][slope_type] for i in range(subj['nbsubj'])]
    else:
        return [subj[i][ch][slope_type][0] for i in range(subj['nbsubj'])]

################################################################################

# Make directory for this run and write parameters to file.
current_time = str(datetime.datetime.now()).split()[0]
export_dir_name = export_dir + current_time + '-' + montage + '/'
num = 1
while os.path.isdir(export_dir_name):
    export_dir_name = export_dir + current_time + '-' + montage + '-' + str(num) + '/'
    num += 1
export_dir = export_dir_name
os.mkdir(export_dir)
params = open(export_dir + 'parameters.txt', 'w')
params.write('Time: ' + str(datetime.datetime.now()))
params.write('\nrecompute_psds = ' + str(recompute_psds))
params.write('\npsd_buffer_lofreq = ' + str(psd_buffer_lofreq))
params.write('\npsd_buffer_hifreq = ' + str(psd_buffer_hifreq))
params.write('\nfitting_func = ' + str(fitting_func))
params.write('\nfitting_lofreq = ' + str(fitting_lofreq))
params.write('\nfitting_hifreq = ' + str(fitting_hifreq))
params.write('\nnwins_upperlimit = ' + str(nwins_upperlimit))
params.write('\nimport_dir = ' + str(import_dir))
params.write('\nexport_dir = ' + str(export_dir))
params.close()

# Compute per-channel PSDs for each subject.
subj = compute_subject_psds(import_dir, '../../../data/auxilliary/ya-oa.csv')
subj['time_computed'] = current_time
filename = export_dir + 'subj-no-fitting.npy'
np.save(filename, subj)

# Select fitting function
if fitting_func == 'linreg':
    regr = linreg_slope
elif fitting_func == 'ransac':
    regr = ransac_slope

# Fit lines to slopes using specified function and frequency range.
print('Fitting to PSD slopes...')
subj = fit_slopes(subj, regr, fitting_lofreq, fitting_hifreq)

# Save results.
filename = export_dir + 'subj-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '-' + fitting_func + '.npy'
subj['time_computed'] = current_time
np.save(filename, subj)

# Define channels, these will form labels for our table.
if montage == 'sensor-level':
    channels = ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30','B31','B32','FRONTAL','LTEMPORAL','CENTRAL','RTEMPORAL','OCCIPITAL']
elif montage == 'dmn':
    channels = ['PCC','PCCr','PCCv','PCCh','mPFC','mPFCr','mPFCv','mPFCh','LAG','LAGr','LAGv','LAGh','RAG','RAGr','RAGv','RAGh','LLatT','LLatTe1','LLatTe2','LLatTe3','RLatT','RLatTe1','RLatTe2','RLatTe3','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise1M','Noise1M1','Noise1M2','Noise1M3','Noise2M','Noise2M1','Noise2M2','Noise2M3']
elif montage == 'frontal':
    channels = ['LdlPFC','LdlPFC1','LdlPFC2','LdlPFC3','RdlPFC','RdlPFC1','RdlPFC2','RdlPFC3','LFRont','LFront1','LFront2','LFront3','RFront','RFront1','RFront2','RFront3','LIPL','LIPLr','LIPLv','LIPLh','RIPL','RIPLr','RIPLv','RIPLh','LIPS','LIPSr','LIPSv','LIPSh','RIPS','RIPSr','RIPSv','RIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','NoiseF','NoiseF1','NoiseF2','NoiseF3']
elif montage == 'dorsal':
    channels = ['LFEF','LFEFr','LFEFv','LFEFh','RFEF','RFEFr','RFEFv','RFEFh','LaIPS','LaIPSr','LaIPSv','LaIPSh','RaIPS','RaIPSr','RaIPSv','RaIPSh','LpIPS','LpIPSr','LpIPSv','LpIPSh','RpIPS','RpIPSr','RpIPSv','RpIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise3L','Noise3L1','Noise3L2','Noise3L3','Noise4R','Noise4R1','Noise4R2','Noise4R3']
elif montage == 'ventral':
    channels = ['LIFG','LIFGr','LIFGv','LIFGh','RIFG','RIFGr','RIFGv','RIFGh','LMFG','LMFGr','LMFGv','LMFGh','RMFG','RMFGr','RMFGv','RMFGh','LTPJ','LTPJr','LTPJv','LTPJh','RTPJ','RTPJr','RTPJv','RTPJh','LSTG','LSTGr','LSTGv','LSTGh','RSTG','RSTGr','RSTGv','RSTGh','NoiseL','NoiseL1','NoiseL2','NoiseL3','NoiseR','NoiseR1','NoiseR2','NoiseR3','NoiseF','NoiseF1','NoiseF2','NoiseF3','Noise','Noise1','Noise2','Noise3']
else:
    raise Exception('ERROR: Montage not recognized.')

data = {}
data['SUBJECT'] = [subj[i]['name'] for i in range(subj['nbsubj'])]
data['CLASS']   = [subj[i]['class'] for i in range(subj['nbsubj'])]
data['AGE']     = [subj[i]['age'] for i in range(subj['nbsubj'])]
data['NWINDOWS_EYESC'] = [subj[i]['eyesc_nwindows'] for i in range(subj['nbsubj'])]
data['NWINDOWS_EYESO'] = [subj[i]['eyeso_nwindows'] for i in range(subj['nbsubj'])]

df = pd.DataFrame(data)
df = df[['SUBJECT', 'CLASS', 'AGE', 'NWINDOWS_EYESC', 'NWINDOWS_EYESO']]

# Add each subject's mean slope.
df['AVG_PSD_EYESC'] = get_subject_slopes(subj, -1, 'eyesC_slope')
df['AVG_PSD_EYESO'] = get_subject_slopes(subj, -1, 'eyesO_slope')

# Now add slopes for every channel from each subject.
for ch in range(len(channels)):
    df[channels[ch] + '_EYESC'] = get_subject_slopes(subj, ch, 'eyesC_slope')
for ch in range(len(channels)):
    df[channels[ch] + '_EYESO'] = get_subject_slopes(subj, ch, 'eyesO_slope')

# Export results
filename = export_dir + 'rs-full-' + montage + '-' + fitting_func + '-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '.csv'
print('Saving fitted slopes at:\n', filename)
df.to_csv(filename, index=False)

