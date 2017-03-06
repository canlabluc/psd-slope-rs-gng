import os
import glob
import seaborn
import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import numpy.fft
import scipy.signal
from scipy.stats import linregress
from sklearn import linear_model


#########################################################
# Data importing functions

def get_filelist(import_path):
    matfiles = []
    for root, dirs, files in os.walk(import_path):
        matfiles += glob.glob(os.path.join(root, '*.mat'))
    return matfiles

def import_subject(subj, i, import_path):
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

#########################################################
# Compute average and per-channel PSDs

def _print_window_info(events, port_code):
    evts = [[events[i][1], events[i+1][1]] for i in range(len(events)) if events[i][0] == port_code]
    total_wins = 0
    total_secs = 0
    for e in evts:
        if (e[1]-e[0]) >= 1024:
            pts  = e[1]-e[0]
            secs = (e[1]-e[0]) // 512
            nwin = (e[1]-e[0]) // 512 - 1
            total_wins += nwin
            total_secs += secs
            print('Event {}:\t{} points, {} seconds, {} windows'.format(e, pts, secs, nwin))
    print('Total windows able to be extracted: ', total_wins)

def get_windows(data, events, port_code, nperwindow=512*2, noverlap=512):
    windows = []
    # The following line restructures events of type port_code into the
    # following format:
    #              [latency, end_of_event]
    evts = [[events[i][1], events[i+1][1]] for i in range(len(events)) if events[i][0] == port_code]
    for event in evts:
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

def compute_subject_psds(import_path):
    """
    Import all subjects and compute per-channel as well as average PSDs.
    """
    matfiles = get_filelist(import_path)

    # Temporary -- I believe these aren't included in the samples-features matrix
    if import_path == '../data/pipeline-full/oaExclFiltCARClust-mat/':
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127132.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127133.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127134.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127140.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127154.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127159.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127160.mat')
        matfiles.remove('../data/pipeline-full/oaExclFiltCARClust-mat/120127167.mat')

    subj = {}
    subj['nbsubj'] = len(matfiles)
    subj['f'] = np.linspace(0, 256, 513)
    subj['f'] = subj['f'].reshape(len(subj['f']), 1)
    subj['f_rm_alpha'] = remove_freq_buffer(subj['f'], 7, 14)
    for i in range(len(matfiles)):
        subj = import_subject(subj, i, matfiles[i])
        for ch in range(subj[i]['nbchan']):
            subj[i][ch] = {}
            eyesC_windows = get_windows(subj[i]['data'][ch], subj[i]['events'], 'C1')
            eyesO_windows = get_windows(subj[i]['data'][ch], subj[i]['events'], 'O1')
            subj[i][ch]['eyesC_psd'] = welch(eyesC_windows, 512)
            subj[i][ch]['eyesO_psd'] = welch(eyesO_windows, 512)
            subj[i][ch]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['eyesC_psd'], 7, 14)
            subj[i][ch]['eyesO_psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['eyesO_psd'], 7, 14)
        subj[i]['data'] = np.nan # No longer needed, so clear it from memory
        subj[i]['eyesC_psd'] = np.mean([subj[i][ch]['eyesC_psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['eyesO_psd'] = np.mean([subj[i][ch]['eyesO_psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i]['eyesC_psd'], 7, 14)
        subj[i]['eyesO_psd_rm_alpha'] = remove_freq_buffer(subj[i]['eyesO_psd'], 7, 14)
        print("Processed: ", subj[i]['name'])
    subj['eyesC_psd'] = np.mean([subj[i]['eyesC_psd'] for i in range(subj['nbsubj'])], axis=0)
    subj['eyesO_psd'] = np.mean([subj[i]['eyesO_psd'] for i in range(subj['nbsubj'])], axis=0)
    return subj

#########################################################
# Computing of per-channel and average PSD slopes

def linreg_slope(f, psd, lofreq, hifreq):
    """
    Fits line to the PSD, using regular linear regression.
    Returns slope and fit line.
    """
    model = linear_model.LinearRegression()
    model.fit(f[2*2:24*2], np.log10(psd[2*2:24*2]))
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
    # Fitting on the grand average PSD of all subjects
    eyesC_slope_and_fitline = regr_func(subj['f'], subj['eyesC_psd'], lofreq, hifreq)
    eyesO_slope_and_fitline = regr_func(subj['f'], subj['eyesO_psd'], lofreq, hifreq)
    subj['eyesC_slope'], subj['eyesC_fitline'] = eyesC_slope_and_fitline
    subj['eyesO_slope'], subj['eyesO_fitline'] = eyesO_slope_and_fitline
    for i in range(subj['nbsubj']):
        # Per-subject PSD average fitting
        eyesC_slope_and_fitline = regr_func(subj['f'], subj[i]['eyesC_psd'], lofreq, hifreq)
        eyesO_slope_and_fitline = regr_func(subj['f'], subj[i]['eyesO_psd'], lofreq, hifreq)
        subj[i]['eyesC_slope'], subj[i]['eyesC_fitline'] = eyesC_slope_and_fitline
        subj[i]['eyesO_slope'], subj[i]['eyesO_fitline'] = eyesO_slope_and_fitline
        for ch in range(subj[i]['nbchan']):
            # Per-channel PSD fitting
            eyesC_slope_and_fitline = regr_func(subj['f'], subj[i][ch]['eyesC_psd_rm_alpha'], lofreq, hifreq)
            eyesO_slope_and_fitline = regr_func(subj['f'], subj[i][ch]['eyesO_psd_rm_alpha'], lofreq, hifreq)
            subj[i][ch]['eyesC_slope'], subj[i][ch]['eyesC_fitline'] = eyesC_slope_and_fitline
            subj[i][ch]['eyesO_slope'], subj[i][ch]['eyesO_fitline'] = eyesO_slope_and_fitline
    return subj
