
# coding: utf-8

# # gonogo pipeline - power estimates across traditional frequency bands
# 
# This notebook covers the power estimates across traditional frequency bands for the Go/NoGo data. Estimates are calculated for both Go and NoGo prompts, and for two different time segments for both prompts:
# - 150-300ms post-stimulus
# - 300-600ms post-stimulus
# 
# Data was pre-processed in MATLAB, using the functions and parameters found in `_psd-slope/data/GNG-power-estimates/pipeline.txt`. The notebook proceeds as follows:
# 
# ##### For both the `GO_PROMPT` and `NOGO_PROMPT` events, we do the following:
# 1. Grab two sets of segments, post-stimulus:
#     **A.** 150ms - 300ms
#     **B.** 300ms - 600ms
# 2. Hamming-window the elements of each set together to construct two continuous recordings, **A** and **B**. 
# 3. Compute the PSD of **A** and **B**.
# 4. Compute power across traditional bands for both. 
# 5. Add information to csv, and write to disk.
# 
# 
# ## NOTES
# 
# At a sampling rate of 512, the highest frequency we can obtain from the PSD is 512/2 = 256 Hz. What happens, however, when we're sampling less than a second?
# 
# We need to discard recordings in which we obtain less than 1024 time points of any of the epochs. We can't really calculate a good PSD from less than this. 

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
from scipy.integrate import simps
from sklearn import linear_model
mpl.rcParams['figure.figsize'] = (16, 10)


# ### Functions

# In[2]:

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
        subj[i]['events'].append([event[0][0].strip(), event[1][0][0], event[2][0][0]])
    subj[i]['data'] = np.squeeze(datafile['data'])
    subj[i]['nbchan'] = len(subj[i]['data'])
    return subj

def get_segments(data, events, port_code, seg_start, nperseg):
    # The following line restructures events of type port_code into the 
    # following format:
    #         [start_time, end_time]
    segments = []
#     print(events)
    evts = [[events[i][1]+seg_start, events[i][1]+seg_start+nperseg] for i in range(len(events)) if events[i][0] == port_code]
    for event in evts:
        segments.append(data[event[0] : event[1]])
    return segments

def welch(segments):
    """ Takes segments grabbed using grab_segments(), calculates each segment's
    PSD, and averages PSDs to return an average PSD.
    """
#     psds = [sp.signal.welch(seg, fs=512, nperseg=len(seg), window='hamming') for seg in segments]
#     print(len(psds))
#     print(len(psds[0]))/
#     return np.mean(psd/s, axis=0)
    # First, hamming-window connect all segments.
    contig = np.concatenate([sp.signal.hamming(len(seg))*seg for seg in segments], axis=0)
#     print(len(contig))
    return sp.signal.welch(contig, 512, window='hamming')[1]

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
    for i in range(len(matfiles)):
        
        subj = import_subject(subj, i, matfiles[i])
        print(subj[i]['name'])
        subj[i]['age']   = df[df.SUBJECT == subj[i]['name']].AGE.values[0]
        subj[i]['class'] = df[df.SUBJECT == subj[i]['name']].CLASS.values[0]
        subj[i]['sex']   = df[df.SUBJECT == subj[i]['name']].SEX.values[0]

        for ch in range(subj[i]['nbchan']):
            subj[i][ch] = {}            
            # Grab this channel's epochA and epochB segments for GO, and compute
            # the mean PSD for the epochA and epochB sets
            segmentsGO_A = get_segments(subj[i]['data'][ch], subj[i]['events'], 'GO_PROMPT', 77, 77)
            segmentsGO_B = get_segments(subj[i]['data'][ch], subj[i]['events'], 'GO_PROMPT', 154, 154)
            subj[i][ch]['psd_GO_A'] = welch(segmentsGO_A)
            subj[i][ch]['psd_GO_B'] = welch(segmentsGO_B)
                        
#             if len(subj[i][ch]['psd_GO_A']) != 513:
#                 print("-------------- {} | {}".format(subj[i]['name'], ch))
#             if len(subj[i][ch]['psd_GO_B']) != 513:
#                 print("-------------- {} | {}".format(subj[i]['name'], ch))
            
#             print("PRINTING")
#             plt.plot(subj['f'], subj[i][ch]['psd_GO_A'])
#             plt.plot(subj['f'], subj[i][ch]['psd_GO_B'])
#             plt.xlim([0, 50])
#             plt.ylim([0, 1])
#             print(subj[i][ch]['psd_GO_A'][10*2])
#             print(subj[i][ch]['psd_GO_A'][20*2])
#             return subj[i][ch]['psd_GO_A'];
            
#             print(len(subj[i][ch]['psd_GO_A']))
#             print(len(subj[i][ch]['psd_GO_B']))
            
            subj[i][ch]['GO_A_DELTA'] = freq_band_power(subj[i][ch]['psd_GO_A'], 0.5, 4)
            subj[i][ch]['GO_A_THETA'] = freq_band_power(subj[i][ch]['psd_GO_A'], 4, 7)
            subj[i][ch]['GO_A_ALPHA'] = freq_band_power(subj[i][ch]['psd_GO_A'], 7, 13)
            subj[i][ch]['GO_A_BETA']  = freq_band_power(subj[i][ch]['psd_GO_A'], 13, 30)
#             subj[i][ch]['GO_A_GAMMA'] = freq_band_power(subj[i][ch]['psd_GO_A'], 30, 45)
            
            subj[i][ch]['GO_B_DELTA'] = freq_band_power(subj[i][ch]['psd_GO_B'], 0.5, 4)
            subj[i][ch]['GO_B_THETA'] = freq_band_power(subj[i][ch]['psd_GO_B'], 4, 7)
            subj[i][ch]['GO_B_ALPHA'] = freq_band_power(subj[i][ch]['psd_GO_B'], 7, 13)
            subj[i][ch]['GO_B_BETA']  = freq_band_power(subj[i][ch]['psd_GO_B'], 13, 30)
#             subj[i][ch]['GO_B_GAMMA'] = freq_band_power(subj[i][ch]['psd_GO_B'], 30, 45)
            
            # Grab this channel's epochA and epochB segments for NOGO, and compute
            # the mean PSD for the epochA and epochB sets
            segmentsNOGO_A = get_segments(subj[i]['data'][ch], subj[i]['events'], 'NOGO_PROMPT', 77, 77)
            segmentsNOGO_B = get_segments(subj[i]['data'][ch], subj[i]['events'], 'NOGO_PROMPT', 154, 154)
#             print(len(segmentsNOGO_A))
#             print(len(segmentsNOGO_A[0]))
            subj[i][ch]['psd_NOGO_A'] = welch(segmentsNOGO_A)
            subj[i][ch]['psd_NOGO_B'] = welch(segmentsNOGO_B)
            
#             print(len(subj[i][ch]['psd_NOGO_A']))
            subj[i][ch]['NOGO_A_DELTA'] = freq_band_power(subj[i][ch]['psd_NOGO_A'], 0.5, 4)
            subj[i][ch]['NOGO_A_THETA'] = freq_band_power(subj[i][ch]['psd_NOGO_A'], 4, 7)
            subj[i][ch]['NOGO_A_ALPHA'] = freq_band_power(subj[i][ch]['psd_NOGO_A'], 7, 13)
            subj[i][ch]['NOGO_A_BETA']  = freq_band_power(subj[i][ch]['psd_NOGO_A'], 13, 30)
#             subj[i][ch]['NOGO_A_GAMMA'] = freq_band_power(subj[i][ch]['psd_NOGO_A'], 30, 45)
            
            subj[i][ch]['NOGO_B_DELTA'] = freq_band_power(subj[i][ch]['psd_NOGO_B'], 0.5, 4)
            subj[i][ch]['NOGO_B_THETA'] = freq_band_power(subj[i][ch]['psd_NOGO_B'], 4, 7)
            subj[i][ch]['NOGO_B_ALPHA'] = freq_band_power(subj[i][ch]['psd_NOGO_B'], 7, 13)
            subj[i][ch]['NOGO_B_BETA']  = freq_band_power(subj[i][ch]['psd_NOGO_B'], 13, 30)
#             subj[i][ch]['NOGO_B_GAMMA'] = freq_band_power(subj[i][ch]['psd_NOGO_B'], 30, 45)

            # Store the number of GO and NOGO prompts that we were able to use
            subj[i]['GO_SEGS'] = len(segmentsGO_A)
            subj[i]['NOGO_SEGS'] = len(segmentsNOGO_A)
            
        subj[i]['data'] = np.nan # No longer needed, so clear it from memory
        # Compute subject's mean PSD for GO epochA, epochB and NOGO epochA, epochB
        subj[i]['GO_A_DELTA'] = np.mean([subj[i][ch]['GO_A_DELTA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_A_THETA'] = np.mean([subj[i][ch]['GO_A_THETA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_A_ALPHA'] = np.mean([subj[i][ch]['GO_A_ALPHA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_A_BETA']  = np.mean([subj[i][ch]['GO_A_BETA']  for ch in range(subj[i]['nbchan'])], axis=0)
#         subj[i]['GO_A_GAMMA'] = np.mean([subj[i][ch]['GO_A_GAMMA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_B_DELTA'] = np.mean([subj[i][ch]['GO_B_DELTA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_B_THETA'] = np.mean([subj[i][ch]['GO_B_THETA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_B_ALPHA'] = np.mean([subj[i][ch]['GO_B_ALPHA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['GO_B_BETA']  = np.mean([subj[i][ch]['GO_B_BETA']  for ch in range(subj[i]['nbchan'])], axis=0)
#         subj[i]['GO_B_GAMMA'] = np.mean([subj[i][ch]['GO_B_GAMMA'] for ch in range(subj[i]['nbchan'])], axis=0)
        
        subj[i]['NOGO_A_DELTA'] = np.mean([subj[i][ch]['NOGO_A_DELTA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_A_THETA'] = np.mean([subj[i][ch]['NOGO_A_THETA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_A_ALPHA'] = np.mean([subj[i][ch]['NOGO_A_ALPHA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_A_BETA']  = np.mean([subj[i][ch]['NOGO_A_BETA']  for ch in range(subj[i]['nbchan'])], axis=0)
#         subj[i]['NOGO_A_GAMMA'] = np.mean([subj[i][ch]['NOGO_A_GAMMA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_B_DELTA'] = np.mean([subj[i][ch]['NOGO_B_DELTA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_B_THETA'] = np.mean([subj[i][ch]['NOGO_B_THETA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_B_ALPHA'] = np.mean([subj[i][ch]['NOGO_B_ALPHA'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['NOGO_B_BETA']  = np.mean([subj[i][ch]['NOGO_B_BETA']  for ch in range(subj[i]['nbchan'])], axis=0)
#         subj[i]['NOGO_B_GAMMA'] = np.mean([subj[i][ch]['NOGO_B_GAMMA'] for ch in range(subj[i]['nbchan'])], axis=0)
        print("Processed: ", subj[i]['name'])
    return subj


# In[3]:

def freq_band_power(psd, lofreq, hifreq):
    """ Returns total power in given frequency range. [uV^2]. 
    Utilizes Simpson's Rule to compute the area under the PSD curve. 
    Note that we do not need to square the result, since sp.signal.welch
    provides us with power, or [uV^2] already.
    """
#     print("Length of psd: {}".format(len(psd)))
#     print("Freqs: {} | {}".format(lofreq, hifreq))
#     print("After:         {}".format(len(psd[round(lofreq*2) : round(lofreq*2)])))
    power = simps(psd[round(lofreq*2) : round(hifreq*2)])
    if power == 0:
        print("FREQ: {} | {}".format(lofreq, hifreq))
        print("PSD: {}".format(psd))
    return simps(psd[round(lofreq*2) : round(hifreq*2)])


# ### Processing

# In[4]:

import_dir_oa = '../../data/GNG-power-estimates/oaExclFiltCARClust-mat/'
import_dir_ya = '../../data/GNG-power-estimates/yaExclFiltCARClust-mat/'
subjoa = compute_subject_psds(import_dir_oa, '../../data/GNG/ya-oa-gng.csv')
subjya = compute_subject_psds(import_dir_ya, '../../data/GNG/ya-oa-gng.csv')


# In[5]:

np.save('../../data/GNG-power-estimates/subjoa.npy', subjoa)
np.save('../../data/GNG-power-estimates/subjya.npy', subjya)


# In[29]:

# Construct matrix
data = {}
oa_names = [subjoa[i]['name'] for i in range(subjoa['nbsubj'])]
ya_names = [subjya[i]['name'] for i in range(subjya['nbsubj'])]
data['SUBJECT'] = np.concatenate([oa_names, ya_names], axis=0)

oa_class = [subjoa[i]['class'] for i in range(subjoa['nbsubj'])]
ya_class = [subjya[i]['class'] for i in range(subjya['nbsubj'])]
data['CLASS'] = np.concatenate([oa_class, ya_class], axis=0)

oa_age = [subjoa[i]['age'] for i in range(subjoa['nbsubj'])]
ya_age = [subjya[i]['age'] for i in range(subjya['nbsubj'])]
data['AGE'] = np.concatenate([oa_age, ya_age], axis=0)

df = pd.DataFrame(data)
df = df[['SUBJECT', 'CLASS', 'AGE']]

channels = ["A1","A2","A3","A4","A5","A6","A7","A8","A10","A11","A12","A13","A14","A15","A16","A17","A18","A21","A22","A23","A24","A25","A26","A27","A29","A30","A31","B1","B2","B3","B4","B5","B6","B8","B9","B10","B11","B12","B13","B14","B17","B18","B19","B20","B21","B22","B23","B24","B26","B27","B28","B29","B30","FRONTAL","LTEMPORAL","CENTRAL","RTEMPORAL","OCCIPITAL"]


# In[30]:

oa = [subjoa[i]['GO_A_DELTA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_A_DELTA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHA_DELTA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_A_THETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_A_THETA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHA_THETA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_A_ALPHA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_A_ALPHA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHA_ALPHA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_A_BETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_A_BETA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHA_BETA'] = np.concatenate([oa, ya], axis=0)
# oa = [subjoa[i]['GO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
# ya = [subjya[i]['GO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
# df['AVG_GO_EPOCHA_GAMMA'] = np.concatenate([oa, ya], axis=0)

oa = [subjoa[i]['GO_B_DELTA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_B_DELTA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHB_DELTA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_B_THETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_B_THETA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHB_THETA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_B_ALPHA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_B_ALPHA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHB_ALPHA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['GO_B_BETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['GO_B_BETA'] for i in range(subjya['nbsubj'])]
df['AVG_GO_EPOCHB_BETA'] = np.concatenate([oa, ya], axis=0)
# oa = [subjoa[i]['GO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
# ya = [subjya[i]['GO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
# df['AVG_GO_EPOCHB_GAMMA'] = np.concatenate([oa, ya], axis=0)

oa = [subjoa[i]['NOGO_A_DELTA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_A_DELTA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHA_DELTA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_A_THETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_A_THETA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHA_THETA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_A_ALPHA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_A_ALPHA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHA_ALPHA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_A_BETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_A_BETA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHA_BETA'] = np.concatenate([oa, ya], axis=0)
# oa = [subjoa[i]['NOGO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
# ya = [subjya[i]['NOGO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
# df['AVG_NOGO_EPOCHA_GAMMA'] = np.concatenate([oa, ya], axis=0)

oa = [subjoa[i]['NOGO_B_DELTA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_B_DELTA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHB_DELTA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_B_THETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_B_THETA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHB_THETA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_B_ALPHA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_B_ALPHA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHB_ALPHA'] = np.concatenate([oa, ya], axis=0)
oa = [subjoa[i]['NOGO_B_BETA'] for i in range(subjoa['nbsubj'])]
ya = [subjya[i]['NOGO_B_BETA'] for i in range(subjya['nbsubj'])]
df['AVG_NOGO_EPOCHB_BETA'] = np.concatenate([oa, ya], axis=0)
# oa = [subjoa[i]['NOGO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
# ya = [subjya[i]['NOGO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
# df['AVG_NOGO_EPOCHB_GAMMA'] = np.concatenate([oa, ya], axis=0)

for ch in range(len(channels)):
    oa = [subjoa[i][ch]['GO_A_DELTA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_A_DELTA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHA_DELTA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_A_THETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_A_THETA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHA_THETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_A_ALPHA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_A_ALPHA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHA_ALPHA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_A_BETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_A_BETA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHA_BETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
#     oa = [subjoa[i][ch]['GO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
#     ya = [subjya[i][ch]['GO_A_GAMMA'] for i in range(subjya['nbsubj'])]
#     df['GO_EPOCHA_GAMMA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_B_DELTA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_B_DELTA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHB_DELTA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_B_THETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_B_THETA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHB_THETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_B_ALPHA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_B_ALPHA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHB_ALPHA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['GO_B_BETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['GO_B_BETA'] for i in range(subjya['nbsubj'])]
    df['GO_EPOCHB_BETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
#     oa = [subjoa[i][ch]['GO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
#     ya = [subjya[i][ch]['GO_B_GAMMA'] for i in range(subjya['nbsubj'])]
#     df['GO_EPOCHB_GAMMA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)

    oa = [subjoa[i][ch]['NOGO_A_DELTA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_A_DELTA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHA_DELTA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_A_THETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_A_THETA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHA_THETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_A_ALPHA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_A_ALPHA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHA_ALPHA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_A_BETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_A_BETA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHA_BETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
#     oa = [subjoa[i][ch]['NOGO_A_GAMMA'] for i in range(subjoa['nbsubj'])]
#     ya = [subjya[i][ch]['NOGO_A_GAMMA'] for i in range(subjya['nbsubj'])]
#     df['NOGO_EPOCHA_GAMMA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_B_DELTA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_B_DELTA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHB_DELTA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_B_THETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_B_THETA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHB_THETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_B_ALPHA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_B_ALPHA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHB_ALPHA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
    oa = [subjoa[i][ch]['NOGO_B_BETA'] for i in range(subjoa['nbsubj'])]
    ya = [subjya[i][ch]['NOGO_B_BETA'] for i in range(subjya['nbsubj'])]
    df['NOGO_EPOCHB_BETA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)
#     oa = [subjoa[i][ch]['NOGO_B_GAMMA'] for i in range(subjoa['nbsubj'])]
#     ya = [subjya[i][ch]['NOGO_B_GAMMA'] for i in range(subjya['nbsubj'])]
#     df['NOGO_EPOCHB_GAMMA_' + channels[ch]] = np.concatenate([oa, ya], axis=0)


# In[31]:

df.to_csv('../../data/GNG-power-estimates/ya-oa-gng-power-estimates-2.csv', index=False)


# In[ ]:

subjoa[0]['psd_NOGO_B']


# In[ ]:




# In[ ]:



