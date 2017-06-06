
# coding: utf-8

# # pipeline - gonogo data, 150-300ms and 300-600ms post-stimulus

# In[1]:

get_ipython().magic('matplotlib inline')
import os
import glob
import datetime
import seaborn as sns
import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import numpy.fft
import scipy.signal
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import linregress
from sklearn import linear_model
mpl.rcParams['figure.figsize'] = (16, 10)


# ## Subject importing & PSD slope calculations

# In[4]:

def get_filelist(import_path):
    matfiles = []
    for root, dirs, files in os.walk(import_path):
        matfiles += glob.glob(os.path.join(root, '*.mat'))
    return matfiles

def import_subject(subj, i, import_path):
    """ 
    Imports a single subject and adds them to the subj
    data structure. Additionally, merges 
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

def _print_window_info(events, port_code):
    evts = [[events[i][1], events[i+1][1]] for i in range(len(events)) if events[i][0] == port_code]
    total_wins = 0
    total_secs = 0
    for e in evts:
        if (e[1] - e[0]) >= 1024:
            pts  = e[1] - e[0]
            secs = (e[1] - e[0])//512
            nwin = (e[1] - e[0])//512 - 1
            total_wins += nwin
            total_secs += secs
            print('Event {}:\t{} points, {} seconds, {} windows'.format(e, pts, secs, nwin))
    print('Total windows able to be extracted: ', total_wins)

def get_windows(data, events, port_code, nperwindow=512*2, noverlap=512):
    windows = []
    # The following line restructures events of type port_code into the 
    # following format:
    #         [start_time, end_time]
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

    subj = {}
    subj['nbsubj'] = len(matfiles)
    subj['f'] = np.linspace(0, 256, 513)
    subj['f'] = subj['f'].reshape(len(subj['f']), 1)
    subj['f_rm_alpha'] = remove_freq_buffer(subj['f'], 7, 14)
    for i in range(len(matfiles)):
        
        subj = import_subject(subj, i, matfiles[i])
        subj[i]['age']   = df[df.SUBJECT == subj[i]['name']].AGE.values[0]
        subj[i]['class'] = df[df.SUBJECT == subj[i]['name']].CLASS.values[0]
        subj[i]['sex']   = df[df.SUBJECT == subj[i]['name']].SEX.values[0]

        for ch in range(subj[i]['nbchan']):
            subj[i][ch] = {}
            windows = get_windows(subj[i]['data'][ch], subj[i]['events'], 'C1')
            subj[i][ch]['psd'] = welch(windows, 512)
            subj[i][ch]['psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['psd'], 7, 14)
        subj[i]['nwindows'] = len(windows)
        subj[i]['data'] = np.nan # No longer needed, so clear it from memory
        subj[i]['psd'] = np.mean([subj[i][ch]['psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['psd_rm_alpha'] = remove_freq_buffer(subj[i]['psd'], 7, 14)
        print("Processed: ", subj[i]['name'])
    subj['psd'] = np.mean([subj[i]['psd'] for i in range(subj['nbsubj'])], axis=0)
    return subj


# In[ ]:



