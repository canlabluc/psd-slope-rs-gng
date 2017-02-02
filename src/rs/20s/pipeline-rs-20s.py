
# coding: utf-8

# # pipeline - 20s eyes closed resting state data, spectral slopes
# 
# This notebook contains the revised analysis for calculating the spectral slope data for resting-state 20s recordings. It seems that in prior analyses, subjects were often inserted into the csv in improper order.

# In[11]:

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
from scipy.stats import linregress, ttest_ind
from sklearn import linear_model
mpl.rcParams['figure.figsize'] = (16, 10)


# ## Parameter selection
# - `recompute_psds`: `True` or `False`, for recomputing subject PSDs or loading previous results.
# - `psd_buffer_lofreq`: Scalar, specifies the lower bound for the PSD buffer that we exclude.
# - `psd_buffer_hifreq`: Scalar, specifies the upper bound for the PSD buffer that we exclude.
# - `fitting_func`: `'linreg'` or `'ransac'`, specifies which function to use for fitting. `'linreg'` is simple linear regression. `'ransac'` is RANSAC, a robust fitting method that ignores outliers.
# - `fitting_lofreq`: Scalar, specifies the lower bound for the PSD fitting range.
# - `fitting_hifreq`: Scalar, specifies the upper bound for the PSD fitting range.
# - `import_dir`: String specifying the directory to import results to.
# - `export_dir`: String specifying the directory to export results to.

# In[2]:

recompute_psds = True
psd_buffer_lofreq = 7
psd_buffer_hifreq = 14
fitting_func = 'ransac'
fitting_lofreq = 2
fitting_hifreq = 24
import_dir = '/Users/jorge/Drive/research/_psd-slope/data/rs-20s/ExclFiltCARClust-mat/'
export_dir = '/Users/jorge/Drive/research/_psd-slope/data/rs-20s/results/'


# ##### Set up workspace, print out parameters to text file...

# In[3]:

current_time = '-'.join('_'.join(str(datetime.datetime.now()).split()).split(':'))[:-7]
export_dir = export_dir + current_time + '/'
os.mkdir(export_dir)
params = open(export_dir + 'parameters.txt', 'w')
params.write('recompute_psds = ' + str(recompute_psds))
params.write('\npsd_buffer_lofreq = ' + str(psd_buffer_lofreq))
params.write('\npsd_buffer_hifreq = ' + str(psd_buffer_hifreq))
params.write('\nfitting_func = ' + str(fitting_func))
params.write('\nfitting_lofreq = ' + str(fitting_lofreq))
params.write('\nfitting_hifreq = ' + str(fitting_hifreq))
params.write('\nexport_dir = ' + str(export_dir))
params.close()


# ## Subject Importing & PSD Calculations
# 
# This section imports subject information and computes PSDs using Welch's method. The algorithm proceeds as follows, for each channel:
# 1. Extract as many clean 2-second eyes closed and eyes open segments from the recording. Segments overlap by 50%.
# 2. Multiply each 2-second segment by a 2-second Hamming window. 
# 3. Compute the discrete Fourier transform of each segment, and average DFT'd segments to arrive at a per-channel PSD.
# 
# The PSD is defined as:
# $$
# PSD = log_{10}(\sum\limits_{n=1}^{N}{N})
# $$
# 
# ##### Function Definitions

# In[4]:

older_adults = ['SA', 'SA_Control', 'MCI', 'MCI_Control']
younger_adults = ['DANE']

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
    subj[i]['data'] = np.squeeze(datafile['data'])
    subj[i]['nbchan'] = len(subj[i]['data'])
    return subj

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
        if subj[i]['class'] in older_adults:
            subj[i]['oa'] = 1
        elif subj[i]['class'] in younger_adults:
            subj[i]['oa'] = 0

        for ch in range(subj[i]['nbchan']):
            subj[i][ch] = {}
            subj[i][ch]['eyesC_psd'] = sp.signal.welch(subj[i]['data'][ch], 512, nperseg=512*2, noverlap=512, window='hamming')[1]
            subj[i][ch]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i][ch]['eyesC_psd'], 7, 14)
        subj[i]['data'] = np.nan # No longer needed, so clear from memory
        subj[i]['eyesC_psd'] = np.mean([subj[i][ch]['eyesC_psd'] for ch in range(subj[i]['nbchan'])], axis=0)
        subj[i]['eyesC_psd_rm_alpha'] = remove_freq_buffer(subj[i]['eyesC_psd'], 7, 14)
        print("Processed: ", subj[i]['name'])
    return subj


# In[5]:

# Import EEG for older and younger adults, compute PSDs
if recompute_psds:
    subj = compute_subject_psds(import_dir, '../../data/RS-20s/ya-oa-20s.csv')
    # Save resulting PSDs
    subj['time_computed'] = current_time
    np.save(export_dir + 'subj-no-fitting.npy', subj); subj = []
else:
    # Use files with pre-computed PSDs and a 7 - 14 Hz buffer
    get_ipython().system('cp /Users/jorge/Drive/research/_psd-slope/data/RS-20s/2016-11-05_17-24-19/subj-no-fitting.npy $export_dir')


# ## Fit to Spectral Slopes
# 
# Now we compute PSD slopes for each channel of each subject, and additionally calculate each subject's mean PSD slope. This is found by fitting to the grand average PSD of each subject.
# 
# ##### Function Definitions

# In[6]:

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
        for ch in range(subj[i]['nbchan']):
            # Per-channel PSD fitting
            subj[i][ch]['eyesC_slope'], subj[i][ch]['eyesC_fitline'] = regr_func(subj['f'], subj[i][ch]['eyesC_psd_rm_alpha'], lofreq, hifreq)
    return subj


# ##### Processing

# In[7]:

# Select fitting function
if fitting_func == 'linreg':
    regr = linreg_slope
elif fitting_func == 'ransac':
    regr = ransac_slope

# Load subject PSDs
subj = np.load(export_dir + '/subj-no-fitting.npy').item()

# Fit lines to slopes using specified function and frequency range
subj = fit_slopes(subj, regr, fitting_lofreq, fitting_hifreq)

# Save results
filename = export_dir + 'subj-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '-' + fitting_func + '.npy'
subj['time_computed'] = current_time
np.save(filename, subj); subj = []


# ## Construct Samples-Features Matrix
# 
# Now we construct the samples-features matrix containing the calculated slopes. We use a table that already contains subject numbers, sex, age, and memory class to start off.
# 
# ##### Function Definitions

# In[8]:

def get_subject_slopes(subj, ch, slope_type):
    """ Returns list of slopes for specified channel of slope_type.
    Arguments:
        subj: The subj data structure.
        ch:   Scalar, channel for which to get list of subject slopes.
        slope_type: String, e.g., 'eyesO_slope' or 'eyesC_slope'
    """
    if ch == -1: # Slope of PSD grand average
        return [subj[i][slope_type]     for i in range(subj['nbsubj'])]
    else:
        return [subj[i][ch][slope_type][0] for i in range(subj['nbsubj'])]


# In[9]:

# Define channels, these will form labels for our table:
channels = ["A01","A02","A03","A04","A05","A06","A07","A08","A09","A10","A11","A12","A13","A14","A15","A16","A17","A18","A19","A20","A21","A22","A23","A24","A25","A26","A27","A28","A29","A30","A31","A32","B01","B02","B03","B04","B05","B06","B07","B08","B09","B10","B11","B12","B13","B14","B15","B16","B17","B18","B19","B20","B21","B22","B23","B24","B25","B26","B27","B28","B29","B30","B31","B32","FRONTAL","LTEMPORAL","CENTRAL","RTEMPORAL","OCCIPITAL"]

# Load subject PSDs with fitted slopes
subj = np.load(filename).item()

# Construct matrix
data = {}
data['SUBJECT'] = [subj[i]['name']  for i in range(subj['nbsubj'])]
data['CLASS']   = [subj[i]['class'] for i in range(subj['nbsubj'])]
data['AGE']     = [subj[i]['age']   for i in range(subj['nbsubj'])]

df = pd.DataFrame(data)
df = df[['SUBJECT', 'CLASS', 'AGE']]

# Add each subject's mean slope.
df['AVG_PSD_EYESC'] = get_subject_slopes(subj, -1, 'eyesC_slope')

# Now add slopes for every channel from each subject.
for ch in range(len(channels)):
    df[channels[ch] + '_EYESC'] = get_subject_slopes(subj, ch, 'eyesC_slope')

# Export results
filename = export_dir + 'ya-oa-rs-20s-' + fitting_func + '-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '-eyesc.csv'
print('Saving fitted slopes at:\n{}'.format(filename))
df.to_csv(filename, index=False); df = []


# ## Compute t-tests, logit, lasso
# 
# And run t-tests, as well as LASSO in order to see if there are any group differences. 
# 
# ### T-Tests
# 
# ##### Younger adults vs older adults

# In[51]:

df = pd.read_csv(filename)
ya = df[df.CLASS.isin(['DANE'])]
oa = df[df.CLASS.isin(['SA_Control', 'MCI_Control'])]

channels = list(df.columns.values)[4:]
for ch in channels:
    result = ttest_ind(ya[ch], oa[ch], equal_var=False)
    if result[1] < 0.05:
        print("{}:\t\t\t{:.2f},\t{:.3f}".format(ch, result.statistic, result.pvalue))


# ##### Older adult controls vs SAs

# In[69]:

oa_control = df[df.CLASS.isin(['SA_Control', 'MCI_Control'])]
sa         = df[df.CLASS.isin(['SA'])]

for ch in channels:
    result = ttest_ind(oa_control[ch], sa[ch], equal_var=False)
    if result[1] < 0.05:
        print("{}:\t\t\t{:.2f},\t{:.3f}".format(ch, result.statistic, result.pvalue))


# ### Logistic Regression, younger adults vs older adult controls

# In[66]:

df_logit = df[df.CLASS.isin(['DANE', 'SA_Control', 'MCI_Control'])]

df_logit['CLASS'] = list(map(lambda x: 0 if x == 'DANE' else 1, df_logit['CLASS']))

cols = list(df_logit.columns.values)
cols.remove('SUBJECT')
cols.remove('CLASS')
cols.remove('AGE')

X = df_logit[cols]
y = df_logit.CLASS


# In[67]:

import warnings                 # sklearn is using a deprecated rand function here,
with warnings.catch_warnings(): # and warnings clutter output
    warnings.simplefilter("ignore")
    resamplings = 2000
    rlogit = linear_model.RandomizedLogisticRegression(n_resampling=resamplings)
    rlogit.fit(X, y)
    print("Features sorted by score, using {} resamplings: ".format(resamplings))
    feature_list = sorted(zip(map(lambda x: round(x, 4), rlogit.scores_), cols), reverse=True)
    for f in feature_list[0:25]: # Adjust this if last feature output is nonzero
        print("{}:\t\t\t{:.2f}".format(f[1], f[0]))


# ### Entire dataset, LASSO for age as interest variable.

# In[68]:

X, y = df[cols], df.AGE

import warnings                 # sklearn is using a deprecated rand function here,
with warnings.catch_warnings(): # and warnings clutter output
    warnings.simplefilter("ignore")
    resamplings = 2000
    rlasso = linear_model.RandomizedLasso(n_resampling=resamplings)
    rlasso.fit(X, y)
    print("Features sorted by score, using {} resamplings: ".format(resamplings))
    feature_list = sorted(zip(map(lambda x: round(x, 4), rlasso.scores_), cols), reverse=True)
    for f in feature_list[0:50]: # Adjust this if last feature output is nonzero
        print("{}:\t\t\t{:.2f}".format(f[1], f[0]))


# In[ ]:



