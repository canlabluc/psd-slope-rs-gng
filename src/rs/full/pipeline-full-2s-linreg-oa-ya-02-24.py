import seaborn
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

def linreg_slope(f, psd, lofreq, hifreq):
    """
    Fits line to the PSD, using regular linear regression.
    Returns slope and fit line.
    """
    model = linear_model.LinearRegression()
    model.fit(f[lofreq*2:hifreq*2], np.log10(psd[lofreq*2:hifreq*2]))
    fit_line = model.predict(f)
    return model.coef_[0] * (10**2), fit_line

def fit_slopes(subj, regr_func, lofreq, hifreq):
    """ 
    Takes subj data structure and fits slopes to each subject's PSDs and mean
    PSD, using regr_func and fitting to datapoints between lofreq and hifreq.
    """
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

def get_subject_slopes(subj, ch, slope_type):
    if ch == -1: # Slope of PSD grand average
        return [subj[i][slope_type]     for i in range(subj['nbsubj'])]
    else:
        return [subj[i][ch][slope_type] for i in range(subj['nbsubj'])]

################################################################################

# Load subject slopes
subjoa = np.load('../data/pipeline-full/subjoa-no-fitting.npy').item()
subjya = np.load('../data/pipeline-full/subjya-no-fitting.npy').item()

# Fit lines to slopes between 30 - 45 Hz using regular linear regression
subjoa = fit_slopes(subjoa, linreg_slope, 2, 24)
subjya = fit_slopes(subjya, linreg_slope, 2, 24)

# Save results
np.save('../data/pipeline-full/subjoa-30-45fit.npy', subjoa)
np.save('../data/pipeline-full/subjya-30-45fit.npy', subjya)

channels = ["A1","A2","A3","A4","A5","A6","A7","A8","A10","A11","A12","A13","A14","A15","A16","A17","A18","A21","A22","A23","A24","A25","A26","A27","A29","A30","A31","B1","B2","B3","B4","B5","B6","B8","B9","B10","B11","B12","B13","B14","B17","B18","B19","B20","B21","B22","B23","B24","B26","B27","B28","B29","B30","FRONTAL","LTEMPORAL","CENTRAL","RTEMPORAL","OCCIPITAL"]

# Read in subject matrix -- just contains subject numbers, age and class.
df = pd.read_csv('../data/ya-oa.csv')
df.SUBJECT = list(map(str, df.SUBJECT))

# Add the overall older adult and younger adult PSD slope.
df['AVG_OA_PSD_EYESC'] = subjoa['eyesC_slope']
df['AVG_OA_PSD_EYESO'] = subjoa['eyesO_slope']
df['AVG_YA_PSD_EYESC'] = subjya['eyesC_slope']
df['AVG_YA_PSD_EYESO'] = subjya['eyesO_slope']

# Add each subject's mean slope. This is found by finding the subject's average PSD and
# fitting a line to it.
oa_avg_psd_eyesc = get_subject_slopes(subjoa, -1, 'eyesC_slope')
ya_avg_psd_eyesc = get_subject_slopes(subjya, -1, 'eyesC_slope')
df['AVG_PSD_EYESC'] = np.concatenate([oa_avg_psd_eyesc, ya_avg_psd_eyesc], axis=0)

oa_avg_psd_eyeso = get_subject_slopes(subjoa, -1, 'eyesO_slope')
ya_avg_psd_eyeso = get_subject_slopes(subjya, -1, 'eyesO_slope')
df['AVG_PSD_EYESO'] = np.concatenate([oa_avg_psd_eyeso, ya_avg_psd_eyeso], axis=0)

# Now add slopes for every channel from each subject.
for ch in range(len(channels)):
    oa_psd_eyesc = get_subject_slopes(subjoa, ch, 'eyesC_slope')
    ya_psd_eyesc = get_subject_slopes(subjya, ch, 'eyesC_slope')
    df[channels[ch] + '_EYESC'] = np.concatenate([oa_psd_eyesc, ya_psd_eyesc], axis=0)
for ch in range(len(channels)):
    oa_psd_eyeso = get_subject_slopes(subjoa, ch, 'eyesO_slope')
    ya_psd_eyeso = get_subject_slopes(subjya, ch, 'eyesO_slope')
    df[channels[ch] + '_EYESO'] = np.concatenate([oa_psd_eyeso, ya_psd_eyeso], axis=0)

df.to_csv('../data/pipeline-full/ya-oa-full-linreg-02-24-eyesc-eyeso-sep.csv', index=False)

