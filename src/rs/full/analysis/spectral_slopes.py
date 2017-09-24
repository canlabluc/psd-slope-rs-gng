#!/usr/local/bin/python3
"""
This script computes the neural noise of each channel in each subject,
and fits a line to the PSD in the specified frequency range. A frequency
range can be specified to ignore when fitting (the 'alpha buffer').

The script works for both sensor-level data and BESA source models.
Change the parameters below appropriately before running. Note that
the script will make a new directory inside of `export_dir` and write
all parameters from that run to parameters.txt.
"""

import os
import sys
import glob
import getopt
import datetime

import numpy as np
import scipy as sp
import pandas as pd
import scipy.io
import scipy.signal
from sklearn import linear_model
from collections import OrderedDict

from subject import Subject

###############################################################################

def get_filelist(import_path, extension):
    """
    Returns list of file paths from import_path with specified extension.
    """
    filelist = []
    for root, dirs, files in os.walk(import_path):
        filelist += glob.glob(os.path.join(root, '*.' + extension))
    return filelist


def get_subject_slopes(subj, ch, slope_type):
    """ Returns list of slopes for specified channel of slope_type.
    Parameters
    ----------
        subj : dict
            Dictionary of Subject objects.
        ch : int
            Channel for which to get list of subject slopes.
        slope_type : str
            Which slope to grab. e.g., 'eyesc_slope' or 'eyeso_slope'
    """
    return [subj[i].psds[ch][slope_type + '_slope'][0] for i in range(subj['nbsubj'])]


def make_export_dir(params):
    """ Creates directory into which we'll save the outputs of the
    analysis.

    Since the same analysis is often run multiple times while
    the script is tweaked, make_export_dir appends the run number at
    the end of the directory name in order to avoid overwriting results
    from a run done previously on the same day. An example of what the
    created directory will look like:
        2017-05-20-sensor-level/
    Parameters
    ----------
        params : dict
            Dictionary containing at least the following keys:
            'export_dir' : str
                Absolute path to the export directory
            'time' : str
                String generated through the following code:
                    str(datetime.datetime.now()).split()[0]
            'montage' : str
                Montage that we're running. e.g., 'sensor-level', or
                'dmn' for default mode network.
    """
    export_dir_name = (
        params['export_dir'] + '/' + params['time'] + '-rs-' + params['montage'] + '/'
    )
    num = 1
    while os.path.isdir(export_dir_name):
        export_dir_name = params['export_dir'] + '/' + params['time'] + '-' +\
                                      params['montage'] + '-' + str(num) + '/'
        num += 1
    params['export_dir'] = export_dir_name
    os.mkdir(params['export_dir'])

    # Write parameters to terminal and to parameters.txt.
    with open(params['export_dir'] + 'parameters.txt', 'w') as params_file:
        print()
        for p in params:
            line = ' {}: {}'.format(p, str(params[p]))
            print(line)
            params_file.write(line + '\n')
        print()
    return params


def get_cmdline_params(params):
    """ Fetches command-line parameters.
    Parameters
    ----------
        params : dict
            Dictionary into which we'll place command-line params.
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hm:i:e:c:o:p:',
            ['fittingfunc=', 'fittinglo=', 'fittinghi=', 'bufferlo=', 'bufferhi=',
             'trialprotocol=', 'nwinsupper='])
    except getopt.GetoptError:
        print('Error: Bad input. To run:\n')
        print('\tspectral_slopes.py -m <montage> -i <import_dir> -o <export_dir> -p <trial_protocol>\n')
        print('Or, manually modify program parameters and run without command-line args.')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Help:\n')
            print('\tspectral_slopes.py -m <montage> -i <import_dir> -o <export_dir> -p <trial_protocol>\n')
            print('Or, manually modify program parameters and run without command-line args.')
            sys.exit(2)

        # Parameters for file imports
        elif opt == '-m':
            params['montage'] = arg
        elif opt == '-i':
            params['import_dir_mat'] = arg
        elif opt == '-e':
            params['import_dir_evt'] = arg
        elif opt == '-c':
            params['import_path_csv'] = arg
        elif opt == '-o':
            params['export_dir'] = arg

        # Parameters for fitting to PSD
        elif opt == '--fittingfunc':
            params['fitting_func'] = arg
        elif opt == '--fittinglo':
            params['fitting_lofreq'] = int(arg)
        elif opt == '--fittinghi':
            params['fitting_hifreq'] = int(arg)
        elif opt == '--bufferlo':
            params['psd_buffer_lofreq'] = int(arg)
        elif opt == '--bufferhi':
            params['psd_buffer_hifreq'] = int(arg)
        elif opt == '--trialprotocol':
            params['trial_protocol'] = arg
        elif opt == '--nwinsupper':
            params['nwins_upperlimit'] = int(arg)
    return params

###############################################################################

def main(argv):
    """
    Parameters : Change these before running, or change them through
    the command-line flags.
    ----------
    montage : string
        Command-line flag: -m
        montage we're running spectral_slopes on. Options are:
            'dmn': Default mode network source model.
            'frontal': Frontal source model.
            'dorsal': Dorsal attention source model.
            'ventral': Ventral attention source model.
            'sensor-level': For running the original sensor-level data.
    psd_buffer_lofreq : float
        Command-line flag: --bufferlo=
        lower frequency bound for the PSD buffer we exclude from fitting.
    psd_buffer_hifreq : float
        Command-line flag: --bufferhi=
        upper frequency bound for the PSD buffer we exclude from fitting.
    fitting_func : string
        Command-line flag: --fittingfunc=
        function we use for fitting to the PSDs. Options are:
            'linreg': Simple linear regression.
            'ransac': RANSAC, a robust fitting method.
    fitting_lofreq : float
        Command-line flag: --fittinglo=
        lower frequency bound for the PSD fitting.
    fitting_hifreq : float
        Command-line flag: --fittinghi=
        higher frequency bound for the PSD fitting.
    trial_protocol : string
        Command-line flag: -p
        specifies whether to modify trial lengths. available options:
            'match_OA': cuts younger adult trials down by half in order
            to make them match older adult trial lengths.
    nwins_upperlimit : int
        upper limit on number of windows to extract from the younger
        adults. A value of 0 means no upper limit.
    import_dir_mat : string
        Command-line flag: -i
        directory from which we import .mat EEG files.
    import_dir_evt : string
        Commane-line flag: -e
        directory from which we import .evt event files.
    export_dir : string
        Command-line flag: -o
        directory to which we export the results, as a .csv file.
    """

    params = OrderedDict()
    params['montage']           = 'sensor-level'
    params['psd_buffer_lofreq'] = 7
    params['psd_buffer_hifreq'] = 14
    params['fitting_func']      = 'ransac'
    params['fitting_lofreq']    = 2
    params['fitting_hifreq']    = 24
    params['trial_protocol']    = 'match_OA'
    params['nwins_upperlimit']  = 0
    params['import_path_csv']   = 'data/auxilliary/ya-oa-have-files-for-all-conds.csv'
    params['import_dir_mat']    = 'data/rs/full/sensor-level/ExclFiltCARClust-mat/'
    params['import_dir_evt']    = 'data/rs/full/evt/clean/'
    params['export_dir']        = 'data/runs/'

    ###########################################################################

    # Generate information about current run.
    params['time'] = str(datetime.datetime.now()).split()[0]
    with open('.git/refs/heads/master', 'r') as f:
        params['commit'] = f.read()[0:7]
    params.move_to_end('commit', last=False)
    params.move_to_end('time', last=False)

    # Take in command-line args, if they are present.
    params = get_cmdline_params(params)

    # Make a directory for this run and write parameters to terminal and file.
    params = make_export_dir(params)

    ##########################################################################
    # Compute PSDs and fit to slopes.

    # Import subject class and age from auxilliary csv.
    matfiles = get_filelist(params['import_dir_mat'], 'mat')
    df = pd.read_csv(params['import_path_csv'])
    df.SUBJECT = df.SUBJECT.astype(str)
    df.CLASS   = df.CLASS.astype(str)
    df.AGE     = df.AGE.astype(int)

    # Check whether we're missing any subject information (i.e., we have the
    # subject EEG, but they're not present in the .csv).
    missing = [f for f in matfiles if f.split('/')[-1][:-4] not in set(df.SUBJECT)]
    for f in missing:
        print('NOTE: Specified csv does not contain information for subject {}'.format(f.split('/')[-1][:-4]))
    if len(missing) > 0:
        print('To include the above subjects, add their information to the .csv file.\n')
    matfiles = [f for f in matfiles if f not in missing]

    # Import EEG data for each subject.
    subj = {}
    subj['nbsubj'] = len(matfiles)
    for i in range(len(matfiles)):

        # Organize subject info into a Subject object.
        subj_name = matfiles[i].split('/')[-1][:-4]
        print('Processing: {}'.format(subj_name))
        group = df[df.SUBJECT == subj_name].CLASS.values[0]
        age   = df[df.SUBJECT == subj_name].AGE.values[0]
        sex   = df[df.SUBJECT == subj_name].SEX.values[0]
        subj[i] = Subject(matfiles[i], params['import_dir_evt'] + subj_name + '.evt', group, age, sex)

        # Modify trial lengths if needed, and compute per-channel PSDs.
        print('Computing PSDs... ', end='')
        if params['trial_protocol'] == 'match_OA' and group == 'DANE':
            subj[i].modify_trial_length(0, 30)
        subj[i].compute_ch_psds(nwins_upperlimit=params['nwins_upperlimit'])
        print('Done.')

        # Fit line to PSD slopes using specified fitting function across
        # specified fitting range with specified exclusion buffer.
        print('Fitting slopes... ', end='')
        if params['fitting_func'] == 'linreg':
            regr = subj[i].linreg_slope
        elif params['fitting_func'] == 'ransac':
            regr = subj[i].ransac_slope
        subj[i].fit_slopes(params['fitting_func'], params['psd_buffer_lofreq'],
                           params['psd_buffer_hifreq'], params['fitting_lofreq'],
                           params['fitting_hifreq'])
        print('Done.\n')

    ##########################################################################
    # Export results to .csv file, save all slopes/PSDs to disk.

    # Construct Pandas dataframe with subject information and slopes.
    data = {}
    data['SUBJECT'] = [subj[i].name for i in range(subj['nbsubj'])]
    data['CLASS']   = [subj[i].group for i in range(subj['nbsubj'])]
    data['AGE']     = [subj[i].age for i in range(subj['nbsubj'])]
    data['NWINDOWS_EYESC'] = [subj[i].nwins_eyesc for i in range(subj['nbsubj'])]
    data['NWINDOWS_EYESO'] = [subj[i].nwins_eyeso for i in range(subj['nbsubj'])]
    df = pd.DataFrame(data)
    df = df[['SUBJECT', 'CLASS', 'AGE', 'NWINDOWS_EYESC', 'NWINDOWS_EYESO']]
    for ch in range(subj[0].nbchan):
        df[subj[0].chans[ch] + '_EYESC'] = get_subject_slopes(subj, ch, 'eyesc')
    for ch in range(subj[0].nbchan):
        df[subj[0].chans[ch] + '_EYESO'] = get_subject_slopes(subj, ch, 'eyeso')

    # Export results to file directory.
    filename = (params['export_dir'] + 'rs-full-' + params['montage'] + '-' +
                params['fitting_func'] + '-' + str(params['fitting_lofreq']) +
                '-' + str(params['fitting_hifreq']) + '.csv')
    print('Saving fitted slopes at:\n', filename)
    df.to_csv(filename, index=False)

    # Write subj dictionary containing fits and slopes to disk. Note
    # that each Subject object is converted to a dictionary for easier
    # future access.
    filename = (params['export_dir'] + 'subj-' + str(params['fitting_lofreq']) +
                '-' + str(params['fitting_hifreq']) + '-' + params['fitting_func'] + '.npy')
    subj['time_computed'] = params['time']
    for i in range(subj['nbsubj']):
        subj[i] = subj[i].__dict__
    np.save(filename, subj)


if __name__ == '__main__':
    main(sys.argv)
