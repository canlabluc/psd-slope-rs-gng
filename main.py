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

from subprocess import check_call

## PARAMETERS #################################################################
###############################################################################

# Analysis and file import parameters
analysis                 = 'rs'
montage                  = 'source-dmn'
import_dir_raw_evt       = 'data/rs/full/evt/raw/'
import_dir_processed_evt = 'data/rs/full/evt/processed/'
import_dir_mat           = 'data/rs/full/source-dmn/MagFiltCAR-mat/'
import_path_csv          = 'data/auxilliary/ya-oa-have-files-for-all-conds.csv'
export_dir               = 'data/runs/'

# PSD Buffer and Fitting Parameters
fitting_func      = 'ransac'
fitting_lofreq    = 2
fitting_hifreq    = 24
psd_buffer_lofreq = 7
psd_buffer_hifreq = 14

# RESTING-STATE-SPECIFIC PARAMETERS
trial_protocol    = 'match_OA'
nwins_upperlimit  = 0

# GO-NOGO-SPECIFIC PARAMETERS
preset_analysis   = 'all_clean_data' # 'fixation_period', 'intertrial_interval', 'custom_marker'
custom_marker     = 'FIXATION'
window_range_lo   = -500
window_range_hi   = 0


## SCRIPT (Do not modify unless you know what you're doing) ###################
###############################################################################

## RESTING STATE ANALYSIS ##
if analysis == 'rs':
  # Pre-process .evts based on analysis user requests
  check_call(['python', 'src/rs/full/preprocessing/cl_evtEMSEPreprocessor.py',
                '-i', import_dir_raw_evt,
                '-o', import_dir_processed_evt])

  # Run spectral_slopes.py
  check_call(['python', 'src/rs/full/analysis/spectral_slopes.py',
    '-m', montage,
    '-i', import_dir_mat,
    '-e', import_dir_processed_evt,
    '-c', import_path_csv,
    '-o', export_dir,
    '--fittingfunc='   + fitting_func,
    '--fittinglo='     + str(fitting_lofreq),
    '--fittinghi='     + str(fitting_hifreq),
    '--bufferlo='      + str(psd_buffer_lofreq),
    '--bufferhi='      + str(psd_buffer_hifreq),
    '--trialprotocol=' + trial_protocol,
    '--nwinsupper='    + str(nwins_upperlimit)])

## GO-NOGO ANALYSIS ##
elif analysis == 'gng':
  # Preprocess .evts based on analysis user requests
  check_call(['python', 'src/gng/preprocessing/cl_evtBESAPreprocessor.py',
                '-i', import_dir_raw_evt,
                '-o', import_dir_processed_evt])

  # Run spectral_slopes.py
  check_call(['python', 'src/gng/analysis/spectral_slopes.py',
    '-m', montage,
    '-i', import_dir_mat,
    '-e', import_dir_processed_evt,
    '-c', import_path_csv,
    '-o', export_dir,
    '--fittingfunc='    + fitting_func,
    '--fittinglo='      + str(fitting_lofreq),
    '--fittinghi='      + str(fitting_hifreq),
    '--bufferlo='       + str(psd_buffer_lofreq),
    '--bufferhi='       + str(psd_buffer_hifreq),
    '--presetanalysis=' + preset_analysis,
    '--custommarker='   + custom_marker,
    '--windowrangelo='  + str(window_range_lo),
    '--windowrangehi='  + str(window_range_hi)])
