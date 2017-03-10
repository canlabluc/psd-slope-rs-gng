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
        subj:       The subj data structure.
        ch:         Scalar, channel for which to get list of subject slopes.
        slope_type: String, e.g., 'eyesO_slope' or 'eyesC_slope'
    """
    return [subj[i].psds[ch][slope_type + '_slope'][0] for i in range(subj['nbsubj'])]

###############################################################################

def main(argv):

    ### Parameters ###
    
    montage = 'dmn'
    recompute_psds = True
    psd_buffer_lofreq = 7
    psd_buffer_hifreq = 14
    fitting_func = 'ransac'
    fitting_lofreq = 2
    fitting_hifreq = 24
    nwins_upperlimit = -1
    cut_recording_length = True
    import_dir_mat = '/Users/jorge/Drive/research/_psd-slope/data/rs/full/source-dmn/MagEvtFiltCAR-mat/'
    import_dir_evt = '/Users/jorge/Drive/research/_psd-slope/data/rs/full/evt/clean/'
    export_dir     = '/Users/jorge/Drive/research/_psd-slope/data/runs/'

    ###
    
    # If present, take in command-line arguments.
    try:
        opts, args = getopt.getopt(argv[1:], 'm:i:o:hc')
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
            montage = arg
        elif opt == '-i':
            import_dir_mat = arg
        elif opt == '-o':
            export_dir = arg
        elif opt == '-c':
            cut_recording_length = True

    # Make directory for this run and write parameters to file.
    current_time = str(datetime.datetime.now()).split()[0]
    export_dir_name = export_dir + '/' + current_time + '-' + montage + '/'
    num = 1
    while os.path.isdir(export_dir_name):
        export_dir_name = export_dir + '/' + current_time + '-' + montage + '-' + str(num) + '/'
        num += 1
    export_dir = export_dir_name
    os.mkdir(export_dir)

    parameters = '''
    Time: {0}
    montage = {1}
    recompute_psds = {2}
    psd_buffer_lofreq = {3}
    psd_buffer_hifreq = {4}
    fitting_func = {5}
    fitting_lofreq = {6}
    fitting_hifreq = {7}
    nwins_upperlimit = {8}
    cut_recording_length = {9}
    import_dir = {10}
    export_dir = {11}
    '''.format(str(datetime.datetime.now()), montage, str(recompute_psds),
               str(psd_buffer_lofreq), str(psd_buffer_hifreq), str(fitting_func),
               str(fitting_lofreq), str(fitting_hifreq), str(nwins_upperlimit),
               str(cut_recording_length), str(import_dir_mat), str(export_dir))
    params = open(export_dir + 'parameters.txt', 'w')
    params.write(parameters)
    params.close()
    print(parameters)

    ### Now for the actual analysis

    # Compute per-channel PSDs for each subject.
    subj = {}
    matfiles = get_filelist(import_dir_mat, 'mat')
    df = pd.read_csv('data/auxilliary/ya-oa.csv')
    df.SUBJECT = df.SUBJECT.astype(str)
    df.CLASS   = df.CLASS.astype(str)
    df.AGE     = df.AGE.astype(int)

    subjects = set(map(lambda x: x.split('/')[-1][:-4], matfiles))
    missing = subjects - set(df.SUBJECT)
    if len(missing) != 0:
        for s in missing:
            print('ERROR: Specified csv does not contain information for subject {}\n'.format(s))
        raise Exception('Missing subject information from csv. Either remove subject file from\n'+
                        'processing pipeline or add subject information to csv file.')

    subj = {}
    subj['nbsubj'] = len(matfiles)
    for i in range(len(matfiles)):
        subj_name = matfiles[i].split('/')[-1][:-4]
        print('Processing: {}'.format(subj_name))

        group = df[df.SUBJECT == subj_name].CLASS.values[0]
        age   = df[df.SUBJECT == subj_name].AGE.values[0]
        sex   = df[df.SUBJECT == subj_name].SEX.values[0]

        subj[i] = Subject(matfiles[i], import_dir_evt + subj_name + '.evt', age, group, sex)
        print('Computing PSDs... ', end='')
        subj[i].compute_ch_psds()
        print('Done.')

        # Select fitting function
        if fitting_func == 'linreg':
            regr = subj[i].linreg_slope
        elif fitting_func == 'ransac':
            regr = subj[i].ransac_slope
        print('Fitting slopes... ', end='')
        subj[i].fit_slopes(fitting_func, 7, 14)
        print('Done.\n')

        filename = export_dir + 'subj-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '-' + fitting_func + '.npy'
        subj['time_computed'] = current_time # TODO
        np.save(filename, subj)

    # Define channels, these will form labels for our table.
    if montage == 'sensor-level':
        channels = ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30','B31','B32','FRONTAL','LTEMPORAL','CENTRAL','RTEMPORAL','OCCIPITAL']
    elif montage == 'dmn':
        channels = ['PCC','PCCr','PCCv','PCCh','mPFC','mPFCr','mPFCv','mPFCh','LAG','LAGr','LAGv','LAGh','RAG','RAGr','RAGv','RAGh','LLatT','LLatTe1','LLatTe2','LLatTe3','RLatT','RLatTe1','RLatTe2','RLatTe3','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise1M','Noise1M1','Noise1M2','Noise1M3','Noise2M','Noise2M1','Noise2M2','Noise2M3']
    elif montage == 'frontal':
        channels = ['LdlPFC','LdlPFC1','LdlPFC2','LdlPFC3','RdlPFC','RdlPFC1','RdlPFC2','RdlPFC3','LFRont','LFront1','LFront2','LFront3','RFront','RFront1','RFront2','RFront3','LIPL','LIPLr','LIPLv','LIPLh','RIPL','RIPLr','RIPLv','RIPLh','LIPS','LIPSr','LIPSv','LIPSh','RIPS','RIPSr','RIPSv','RIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','NoiseF','NoiseF1','NoiseF2','NoiseF3']
    elif montage == 'dorsal':
        channels = ['LFEF','LFEFr','LFEFv','LFEFh','RFEF','RFEFr','RFEFv','RFEFh','LaIPS','LaIPSr','LaIPSv','LaIPSh','RaIPS','RaIPSr','RaIPSv','RaIPSh','LpIPS','LpIPSr','LpIPSv','LpIPSh','RpIPS','RpIPSr','RpIPSv','RpIPSh','Noise1L','Noise1L1','Noise1L2','Noise1L3','Noise1R','Noise1R1','Noise1R2','Noise1R3','Noise2L','Noise2L1','Noise2L2','Noise2L3','Noise2R','Noise2R1','Noise2R2','Noise2R3','Noise3L','Noise3L1','Noise3L2','Noise3L3','Noise4R','Noise4R1','Noise4R2','Noise4R3']
    elif montage == 'ventral':
        channels = ['LIFG','LIFGr','LIFGv','LIFGh','RIFG','RIFGr','RIFGv','RIFGh','LMFG','LMFGr','LMFGv','LMFGh','RMFG','RMFGr','RMFGv','RMFGh','LTPJ','LTPJr','LTPJv','LTPJh','RTPJ','RTPJr','RTPJv','RTPJh','LSTG','LSTGr','LSTGv','LSTGh','RSTG','RSTGr','RSTGv','RSTGh','NoiseL','NoiseL1','NoiseL2','NoiseL3','NoiseR','NoiseR1','NoiseR2','NoiseR3','NoiseF','NoiseF1','NoiseF2','NoiseF3','Noise','Noise1','Noise2','Noise3']
    else:
        raise Exception('ERROR: Montage not recognized.')

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

    # Export results
    filename = export_dir + 'rs-full-' + montage + '-' + fitting_func + '-' + str(fitting_lofreq) + '-' + str(fitting_hifreq) + '.csv'
    print('Saving fitted slopes at:\n', filename)
    df.to_csv(filename, index=False)


if __name__ == '__main__':
    main(sys.argv[:])
