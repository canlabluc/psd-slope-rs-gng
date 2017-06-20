#!/usr/local/bin/python3
"""
Computes PSDs and fits line to specified frequency range in PSD for
the sensor-level Go-NoGo data.
Change the parameters below prior to running.
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
    Arguments:
        subj:       Dictionary of Subject objects.
        ch:         Scalar, channel for which to get list of subject slopes.
        slope_type: String, e.g., 'eyesc_slope' or 'eyeso_slope'
    """
    return [subj[i].psds[ch][slope_type + '_slope'][0] for i in range(subj['nbsubj'])]

###############################################################################

def main(argv):
    """
    Parameters : Change these before running.
    ----------
    montage : string
        montage we're running spectral_slopes on. Options are:
            'dmn': Default mode network source model.
            'frontal': Frontal source model.
            'dorsal': Dorsal attention source model.
            'ventral': Ventral attention source model.
            'sensor-level': For running the original sensor-level data.

    psd_buffer_lofreq : float
        lower frequency bound for the PSD buffer we exclude from fitting.

    psd_buffer_hifreq : float
        upper frequency bound for the PSD buffer we exclude from fitting.

    fitting_func : string
        function we use for fitting to the PSDs. Options are:
            'linreg': Simple linear regression.
            'ransac': RANSAC, a robust fitting method.

    fitting_lofreq : float
        upper frequency bound for the PSD fitting.

    fitting_hifreq : float
        lower frequency bound for the PSD fitting.

    nwins_upperlimit : int
        upper limit on number of windows to extract from the younger
        adults. A value of 0 means no upper limit.

    protocol : string
        specifies
        specifies whether to cut younger adult trials down by half in
        order to make them match older adult trial lengths.

    import_dir_mat : string
        directory from which we import .mat files.

    import_dir_evt : string
        directory from which we import .evt files.

    export_dir : string
        directory to which we export the results .csv file.
    """

    params = OrderedDict()
    params['recompute_psds']    = True
    params['psd_buffer_lofreq'] = 7
    params['psd_buffer_hifreq'] = 14
    params['fitting_func']      = 'ransac'
    params['fitting_lofreq'] = 14
    params['fitting_hifreq'] = 34
    params['nwins_upperlimit'] = 0
    params['protocol'] = 'all_clean' # 'all_clean_space', 'correct_go_trials', 'correct_nogo_trials',
    params['import_path_csv'] = 'data/auxilliary/ya-oa.csv'
    params['import_dir_mat'] = 'data/gng/ExclFiltCARClust-mat/'
    params['import_dir_evt'] = 'data/gng/evt/clean/'
    params['export_dir']     = 'data/runs/'

    ###########################################################################

    # Make sure we're working at the project root.
    project_path = os.getcwd()
    os.chdir(project_path[:project_path.find('psd-slope') + len('psd-slope')] + '/')

    # Generate information about current run.
    params['Time'] = str(datetime.datetime.now()).split()[0]
    with open('.git/refs/heads/master', 'r') as f:
        params['commit'] = f.read()[0:7]
    params.move_to_end('commit', last=False)
    params.move_to_end('Time', last=False)

    # Take in command-line args, if they are present.
    try:
        opts, args = getopt.getopt(argv[1:], 'i:e:o:h')
    except getopt.GetoptError:
        print('Error: Bad input. To run:\n')
        print('\tgng_spectral_slopes.py -m <montage> -i <import_dir> -o <export_dir>\n')
        print('Or, manually modify program parameters and run without command-line args.')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Help:\n')
            print('\tgng_spectral_slopes.py -m <montage> -i <import_dir> -o <export_dir>\n')
            print('Or, manually modify program parameters and run without command-line args.')
            sys.exit(2)
        elif opt == '-i':
            params['import_dir_mat'] = arg
        elif opt == '-e':
            params['import_dir_evt'] = arg
        elif opt == '-o':
            params['export_dir'] = arg

    # Make a directory for this run and write parameters to terminal and file.
    export_dir_name = params['export_dir'] + '/' + params['Time'] + '/'
    num = 1
    while os.path.isdir(export_dir_name):
        export_dir_name = params['export_dir'] + '/' + params['Time'] + '-' + str(num) + '/'
        num += 1
    params['export_dir'] = export_dir_name
    os.mkdir(params['export_dir'])

    with open(params['export_dir'] + 'parameters.txt', 'w') as params_file:
        print()
        for p in params:
            line = ' {}: {}'.format(p, str(params[p]))
            print(line)
            params_file.write(line + '\n')
        print()

    ##########################################################################
    # Compute PSDs and fit to slopes.

    subj = {}
    matfiles = get_filelist(params['import_dir_mat'], 'mat')
    df = pd.read_csv(params['import_path_csv'])
    df.SUBJECT = df.SUBJECT.astype(str)
    df.CLASS   = df.CLASS.astype(str)
    df.AGE     = df.AGE.astype(int)

    subjects = set(map(lambda x: x.split('/')[-1][:-4], matfiles))
    missing = subjects - set(df.SUBJECT)
    if len(missing) != 0:
        for s in missing:
            print('ERROR: Specified csv does not contain information for subject {}'.format(s))
        raise Exception('\nMissing subject information from csv. Either remove subject file from\n'+
                        'processing pipeline or add subject information to csv file.')

    subj = {}
    subj['nbsubj'] = len(matfiles)
    for i in range(len(matfiles)):

        #######################################
        # Import the subject
        subj_name = matfiles[i].split('/')[-1][:-4]
        print('Processing: {}'.format(subj_name))
        group = df[df.SUBJECT == subj_name].CLASS.values[0]
        age   = df[df.SUBJECT == subj_name].AGE.values[0]
        sex   = df[df.SUBJECT == subj_name].SEX.values[0]
        subj[i] = Subject(matfiles[i], params['import_dir_evt'] + subj_name + '.evt', group, age, sex)

        #######################################
        # Compute PSDs

        if params['protocol'] == 'correct_go_trials':






        print('Computing PSDs... ', end='')
        # TODO
        if params['protocol'] == 'correct_go_trials':
            subj[i].modify_trial_heirarchy() # Event-file reordering should occur inside of the Subject object.
        elif params['protocol'] == 'correct_nogo_trials':
            subj[i].modify_trial_heirarchy()
        subj[i].compute_ch_psds(nwins=params['nwins_upperlimit'])
        print('Done.')

        #######################################
        # Fit to PSD slopes
        print('Fitting slopes... ', end='')
        if params['fitting_func'] == 'linreg':
            regr = subj[i].linreg_slope
        elif params['fitting_func'] == 'ransac':
            regr = subj[i].ransac_slope
        subj[i].fit_slopes(params['fitting_func'], params['psd_buffer_lofreq'],
                           params['psd_buffer_hifreq'], params['fitting_lofreq'],
                           params['fitting_hifreq'])
        print('Done.\n')

    filename = (params['export_dir'] + 'subj-' + str(params['fitting_lofreq']) +
                '-' + str(params['fitting_hifreq']) + '-' + params['fitting_func'] + '.npy')
    subj['time_computed'] = params['Time']
    np.save(filename, subj)

    ##########################################################################
    # Construct Pandas dataframe and export results to .csv file.

    # Construct Pandas dataframe with subject information and slopes.
    data = {}
    data['SUBJECT'] = [subj[i].name for i in range(subj['nbsubj'])]
    data['CLASS']   = [subj[i].group for i in range(subj['nbsubj'])]
    data['AGE']     = [subj[i].age for i in range(subj['nbsubj'])]
    data['NWINDOWS_EYESC'] = [subj[i].nwins_eyesc for i in range(subj['nbsubj'])]
    data['NWINDOWS_EYESO'] = [subj[i].nwins_eyeso for i in range(subj['nbsubj'])]
    df = pd.DataFrame(data)
    df = df[['SUBJECT', 'CLASS', 'AGE', 'NWINDOWS_EYESC', 'NWINDOWS_EYESO']]
    for ch in range(len(channels)):
        df[channels[ch] + '_EYESC'] = get_subject_slopes(subj, ch, 'eyesc')
    for ch in range(len(channels)):
        df[channels[ch] + '_EYESO'] = get_subject_slopes(subj, ch, 'eyeso')

    # Export results to file directory.
    filename = (params['export_dir'] + 'rs-full-' + '-' +
                params['fitting_func'] + '-' + str(params['fitting_lofreq']) +
                '-' + str(params['fitting_hifreq']) + '.csv')
    print('Saving fitted slopes at:\n', filename)
    df.to_csv(filename, index=False)


if __name__ == '__main__':
    main(sys.argv[:])
