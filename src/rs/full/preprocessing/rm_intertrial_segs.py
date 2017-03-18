"""
Final step in processing. Trims clean segments if they partially sit in the
intertrial space, or outright removes them if they wholly sit in intertrial
space.
"""

import os
import sys
import glob

import numpy as np
import pandas as pd


def get_filelist(import_path, extension):
    """
    Returns list of file paths from import_path with specified extension.
    """
    filelist = []
    for root, dirs, files in os.walk(import_path):
        filelist += glob.glob(os.path.join(root, '*.' + extension))
        return filelist


def get_previous_trial(code_latency, trials):
    """
    Returns closest previous trial to the provided latency from the dataframe
    trials.
    """
    prev_trial = trials[trials.Latency < code_latency].tail(1)
    return prev_trial.set_index(np.arange(0, prev_trial.shape[0], 1))


def get_next_trial(code_latency, trials):
    """
    Returns the closest next trial to the provided latency from the dataframe
    trials.
    """
    next_trial = trials[trials.Latency > code_latency].head(1)
    return next_trial.set_index(np.arange(0, next_trial.shape[0], 1))


def in_intertrial(code_latency, trials):
    """
    Returns boolean specifying whether the given latency sits inside of the
    intertrial space defined by the dataframe trials.
    """
    try:
        if get_previous_trial(code_latency, trials).Type.all()[0:2] in TRIALS_STOP:
            return True
        else:
            return False
    except:
        return True


def print_seg_code_information(df, i, trials, error_type):
    try:
        ptrial_code = get_previous_trial(df.iloc[i].Latency, trials).Type.all()
        ptrial_latency = get_previous_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
    except:
        ptrial_code = 'SOF'
        ptrial_latency = 0
    try:
        ntrial_code = get_next_trial(df.iloc[i].Latency, trials).Type.all()
        ntrial_latency = get_next_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
    except:
        ntrial_code = 'EOF'
        ntrial_latency = -1
    print(error_type + ': {} | PREV: {}:{:0.2f} | {}:{:0.2f}->{}:{:0.2f} | NEXT: {}:{:0.2f}'.format(
        fname,
        ptrial_code,
        ptrial_latency/512,
        df.iloc[i].Type,
        df.iloc[i].Latency/512,
        df.iloc[i+1].Type,
        df.iloc[i+1].Latency/512,
        ntrial_code,
        ntrial_latency/512
    ))


importpath = sys.argv[1]
exportpath = sys.argv[2]

SEGS_START = ['C1', 'O1']
SEGS_STOP = ['C2', 'O2']
TRIALS_START = ['11', '10']
TRIALS_STOP = ['21', '20']
files = get_filelist(importpath, 'evt')
for f in files:
    fname = f.split('/')[-1]
    df = pd.read_csv(f, sep='\t')
    df.Type = df.Type.astype(str)
    df.Latency = df.Latency.astype(int)

    trials = df[[code[0:2] in TRIALS_START + TRIALS_STOP for code in df.Type]]
    found_bad_segs = False

    # CASE 1: Remove segments that fall completely in intertrial space.
    bad_idx = []
    for i in range(df.shape[0]):
        if df.iloc[i].Type in SEGS_START and (in_intertrial(df.iloc[i].Latency, trials) and
                                              in_intertrial(df.iloc[i+1].Latency, trials)):
            print_seg_code_information(df, i, trials, 'INTERTRIAL')
            bad_idx.append(i)
            bad_idx.append(i+1)
            found_bad_segs = True
    df = df.drop(df.index[bad_idx])
    df = df.set_index(np.arange(0, df.shape[0], 1))

    # Trim clean segments that sit partially in intertrial space.
    for i in range(df.shape[0]):

        # CASE 2: Segment head sits inside of the intertrial space. Trim the front of the segment.
        if df.iloc[i].Type[0:2] in SEGS_START and in_intertrial(df.iloc[i].Latency, trials):
            print_seg_code_information(df, i, trials, 'TRIM HEAD')
            df.at[i, 'Latency'] = get_next_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
            found_bad_segs = True

        # CASE 3: Segment tail sits inside of the intertrial space. Trim the end of the segment.
        if df.iloc[i].Type[0:2] in SEGS_STOP and in_intertrial(df.iloc[i].Latency, trials):
            print_seg_code_information(df, i, trials, 'TRIM TAIL')
            df.at[i, 'Latency'] = get_previous_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
            found_bad_segs = True

    if found_bad_segs:
        print()

    df.to_csv(exportpath + fname, sep='\t', index=False)
