""" main.py
Top-level script for running the spectral slopes analysis on the GoNoGo or
resting state data.

To run an analysis, modify the parameters below to the appropriate
values specific to your run. Note that main.ipynb contains the same
code, but in the form of a Jupyter notebook.

PARAMETER OPTIONS
-----------------

analysis : str
  Specifies which analysis to run. Options: 'rs' and 'gng'

preprocess_evts : bool
  Whether to preprocess evts or not. If .evt files have already been
  preprocessed, then this step is not necessary. Note that setting
  this to true causes the script to prompt the user for a path to the
  raw evt files.

import_dir_evt : str
  Path to directory containing the preprocessed .evt files. Path can be
  relative or absolute.

import_dir_mat : str
  Path to directory containing the preprocessed .mat files. These are
  obtained through canlab's EEGLAB plugin.

import_path_csv : str
  Path to csv file containing subjects that we'll be running. The csv
  should have four columns:
    SUBJECT    CLASS    AGE    SEX
  NOTE: The script will ignore subject .mat and .evt files which do not
  have a matching entry in the .csv file.

export_dir : str
  Path to directory in which the script will output results. Note that
  the script will create another directory in export_dir in the format:
    <DATE>-<ANALYSIS>-<MONTAGE>
  Example:
    2017-01-25-gng-sensor-level
  If a directory with that name already exists in export_dir,
  spectral_slopes.py will create another directory with the same name
  and append a number to it, e.g. 2017-01-25-gng-sensor-level-1.

fitting_func : str
  Function the script will use to fit a line to the PSDs. Current
  options: 'ransac', 'linreg' for RANSAC and linear regression,
  respectively.

fitting_lofreq : int
  Lower frequency bound for the PSD fitting.

fitting_hifreq : int
  Upper frequency bound for the PSD fitting.

psd_buffer_lofreq : int
  Lower frequency bound for the PSD exclusion buffer.

psd_buffer_hifreq : int
  Upper frequency bound for the PSD exlcusion buffer.

trial_protocol : str
  Specifies whether or not to modify trial lengths. Available options:
    'match_OA': cuts younger adult trials down by half in order to make
    them match older adult trial lengths.

nwins_upperlimit : int
  Upper limit on number of windows to extract from the younger
  adults. A value of 0 means no upper limit.
"""

import os
import sys
import getopt
from subprocess import check_call

## PARAMETERS #################################################################

# Analysis and file import parameters
analysis        = 'rs'
montage         = 'sensor-level'
preprocess_evts = True
raw_evt_dir     = 'data/rs/full/evt/raw/' # Only used if preprocess_evts is True
import_dir_evt  = 'data/rs/full/evt/clean/'
import_dir_mat  = 'data/rs/full/sensor-level/ExclFiltCARClust-mat/'
import_path_csv = 'data/auxilliary/ya-oa-have-files-for-all-conds.csv'
export_dir      = 'data/runs/'

# Fitting/PSD parameters
fitting_func      = 'ransac'
fitting_lofreq    = 2
fitting_hifreq    = 24
psd_buffer_lofreq = 7
psd_buffer_hifreq = 14
trial_protocol    = 'match_OA'
nwins_upperlimit  = 0

## SCRIPT (Do not modify unless you know what you're doing) ##################
##############################################################################

# RESTING STATE ANALYSIS
if analysis == 'rs':
  # Preprocess .evt files, if they need to be preprocessed.
  if preprocess_evts:
    import_dir_evt_raw = 'data/rs/full/evt/raw/'
    check_call(['python', 'src/rs/full/preprocessing/evt_preprocessing.py',
                  '-i', raw_evt_dir,
                  '-o', import_dir_evt])

  # Run spectral_slopes.py
  check_call(['python', 'src/rs/full/analysis/spectral_slopes.py',
    '-m', montage,
    '-i', import_dir_mat,
    '-e', import_dir_evt,
    '-c', import_path_csv,
    '-o', export_dir,
    '--fittingfunc='   + fitting_func,
    '--fittinglo='     + str(fitting_lofreq),
    '--fittinghi='     + str(fitting_hifreq),
    '--bufferlo='      + str(psd_buffer_lofreq),
    '--bufferhi='      + str(psd_buffer_hifreq),
    '--trialprotocol=' + trial_protocol,
    '--nwinsupper='    + str(nwins_upperlimit)])

# GO-NOGO ANALYSIS
elif analysis == 'gng':
  print('GNG')
