#!/usr/local/bin/python3
"""
Computes PSDs and fits line to specified frequency range in PSD for
BESA source models.
Change parameters below in order to run different datasets.
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
    Parameters
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
    match_OA_protocol : bool
        specifies whether to cut younger adult trials down by half in
        order to make them match older adult trial lengths.
    nwins_upperlimit : int
        upper limit on number of windows to extract from the younger
        adults. A value of 0 means no upper limit.
    import_dir_mat : string
        directory from which we import .mat files.
    import_dir_evt : string
        directory from which we import .evt files.
    export_dir : string
        directory to which we export the results .csv file.
    """

    params = OrderedDict()
    params['montage'] = 'dmn'
    params['recompute_psds'] = True
    params['psd_buffer_lofreq'] = 7
    params['psd_buffer_hifreq'] = 14
    params['fitting_func'] = 'ransac'
    params['fitting_lofreq'] = 14
    params['fitting_hifreq'] = 34
    params['trial_protocol'] = 'match_OA'
    params['nwins_upperlimit'] = 0
    params['import_dir_mat'] = 'data/rs/full/source-dmn/MagCleanEvtFiltCAR-mat/'
    params['import_dir_evt'] = 'data/rs/full/evt/clean/'
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
        opts, args = getopt.getopt(argv[1:], 'm:i:o:hp:')
    except getopt.GetoptError:
        print('Error: Bad input. To run:\n')
        print('\tspectral_slopes.py -m <montage> -i <import_dir> -o <export_dir>\n')
        print('Or, manually modify program parameters and run without command-line args.')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Help:\n')
            print('\tspectral_slopes.py -m <montage> -i <import_dir> -o <export_dir>\n')
            print('Or, manually modify program parameters and run without command-line args.')
            sys.exit(2)
        elif opt == '-m':
            params['montage'] = arg
        elif opt == '-i':
            params['import_dir_mat'] = arg
        elif opt == '-o':
            params['export_dir'] = arg
        elif opt == '-p':
            params['trial_protocol'] = arg

    # Make a directory for this run and write parameters to terminal and file.
    export_dir_name = params['export_dir'] + '/' + params['Time'] + '-' +\
                                                     params['montage'] + '/'
    num = 1
    while os.path.isdir(export_dir_name):
        export_dir_name = params['export_dir'] + '/' + params['Time'] + '-' +\
                                      params['montage'] + '-' + str(num) + '/'
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
    df = pd.read_csv('data/auxilliary/ya-oa.csv')
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

        subj_name = matfiles[i].split('/')[-1][:-4]
        print('Processing: {}'.format(subj_name))
        group = df[df.SUBJECT == subj_name].CLASS.values[0]
        age   = df[df.SUBJECT == subj_name].AGE.values[0]
        sex   = df[df.SUBJECT == subj_name].SEX.values[0]
        subj[i] = Subject(matfiles[i], params['import_dir_evt'] + subj_name + '.evt', group, age, sex)

        print('Computing PSDs... ', end='')
        if params['trial_protocol'] == 'match_OA' and group == 'DANE':
            subj[i].modify_trial_length(0, 30)
        subj[i].compute_ch_psds(nwins_upperlimit=params['nwins_upperlimit'])
        print('Done.')

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

    # Define channel labels for the montage.
    if params['montage'] == 'sensor-level':
        channels = ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30','B31','B32','FRONTAL','LTEMPORAL','CENTRAL','RTEMPORAL','OCCIPITAL']
    elif params['montage'] == 'dmn':
        channels = ['PCC','PCCr','PCCv','PCCh','mPFC','mPFCr','mPFCv','mPFCh','LAG','LAGr','LAGv','LAGh','RAG','RAGr','RAGv','RAGh','LLatT','LLatTe1','LLatTe2','LLatTe3','RLatT','RLatTe1','RLatTe2','RLatTe3','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise1M','Noise1M1','Noise1M2','Noise1M3','Noise2M','Noise2M1','Noise2M2','Noise2M3']
    elif params['montage'] == 'frontal':
        channels = ['LdlPFC','LdlPFC1','LdlPFC2','LdlPFC3','RdlPFC','RdlPFC1','RdlPFC2','RdlPFC3','LFRont','LFront1','LFront2','LFront3','RFront','RFront1','RFront2','RFront3','LIPL','LIPLr','LIPLv','LIPLh','RIPL','RIPLr','RIPLv','RIPLh','LIPS','LIPSr','LIPSv','LIPSh','RIPS','RIPSr','RIPSv','RIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','NoiseF','NoiseF1','NoiseF2','NoiseF3']
    elif params['montage'] == 'dorsal':
        channels = ['LFEF','LFEFr','LFEFv','LFEFh','RFEF','RFEFr','RFEFv','RFEFh','LaIPS','LaIPSr','LaIPSv','LaIPSh','RaIPS','RaIPSr','RaIPSv','RaIPSh','LpIPS','LpIPSr','LpIPSv','LpIPSh','RpIPS','RpIPSr','RpIPSv','RpIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise3L','Noise3L1','Noise3L2','Noise3L3','Noise4R','Noise4R1','Noise4R2','Noise4R3']
    elif params['montage'] == 'ventral':
        channels = ['LIFG','LIFGr','LIFGv','LIFGh','RIFG','RIFGr','RIFGv','RIFGh','LMFG','LMFGr','LMFGv','LMFGh','RMFG','RMFGr','RMFGv','RMFGh','LTPJ','LTPJr','LTPJv','LTPJh','RTPJ','RTPJr','RTPJv','RTPJh','LSTG','LSTGr','LSTGv','LSTGh','RSTG','RSTGr','RSTGv','RSTGh','NoiseL','NoiseL1','NoiseL2','NoiseL3','NoiseR','NoiseR1','NoiseR2','NoiseR3','NoiseF','NoiseF1','NoiseF2','NoiseF3','Noise','Noise1','Noise2','Noise3']
    else:
        raise Exception('ERROR: Montage not recognized.')

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
    filename = (params['export_dir'] + 'rs-full-' + params['montage'] + '-' +
                params['fitting_func'] + '-' + str(params['fitting_lofreq']) +
                '-' + str(params['fitting_hifreq']) + '.csv')
    print('Saving fitted slopes at:\n', filename)
    df.to_csv(filename, index=False)


if __name__ == '__main__':
    main(sys.argv[:])
