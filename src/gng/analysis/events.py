#!/usr/local/bin/python3
"""
Handles event information for gng_spectral_slopes.py.
"""
import numpy as np
import pandas as pd


def get_all_clean_space(df):
    df = df[df.Trigger.isin(['BLINK1', 'BLINK2', 'ARTFCT1', 'ARTFCT2'])]
    df = df.set_index(np.arange(0, df.shape[0], 1))

    # We're going to construct a dataframe that only contains event codes
    # marking clean segments of the recording. We start with creating an empty
    # dataframe with columns for Latency and Trigger.
    clean = pd.DataFrame(columns=['Latency', 'Trigger'])

    # Drop all recording data prior to the beginning of trials.
    clean_idx = 0
    for j in range(df.shape[0]):
        if df.iloc[j].Trigger == 'ARTFCT2':
            clean_idx = j
            break
    clean = clean.append({'Latency': df.iloc[clean_idx].Latency,
                          'Trigger': 'C1'}, ignore_index=True)
    df = df.drop(np.arange(0, clean_idx+1, 1))

    # Grab all segments of recording that occur between blinks and artifacts.
    # These segments are added to the clean dataframe.
    for j in range(df.shape[0]):
        if df.iloc[j].Trigger == 'BLINK1':
            clean = clean.append({'Latency': df.iloc[j].Latency,
                                  'Trigger': 'C2'}, ignore_index=True)
            clean = clean.append({'Latency': df.iloc[j+1].Latency,
                                  'Trigger': 'C1'}, ignore_index=True)
        if df.iloc[j].Trigger == 'ARTFCT1':
            clean = clean.append({'Latency': df.iloc[j].Latency,
                                  'Trigger': 'C2'}, ignore_index=True)
            for k in range(j, df.shape[0]):
                if df.iloc[k].Trigger == 'ARTFCT2':
                    clean = clean.append({'Latency': df.iloc[k].Latency,
                                          'Trigger': 'C1'}, ignore_index=True)
                    break

    # Of the clean segments that were grabbed, drop all which are shorter than
    # two seconds.
    bad_idx = []
    for j in range(clean.shape[0] - 1):
        if (clean.iloc[j].Trigger == 'C1' and clean.iloc[j+1].Trigger == 'C2') and\
           (clean.iloc[j+1].Latency - clean.iloc[j].Latency < 1024):
           bad_idx.append(j); bad_idx.append(j+1)
    clean = clean.drop(bad_idx)
    clean = clean.set_index(np.arange(0, clean.shape[0], 1))

    if clean.iloc[-1].Trigger == 'C1':
        clean = clean.drop(clean.shape[0] - 1)

    subj_name = evt_files[i].split('/')[-1][:-4]
    return clean


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
        if get_previous_trial(code_latency, trials).Type.all() in ['12', '02']:
            return True
        else:
            return False
    except:
        return True


def rm_intertrial_segs(df):
    """
    Checks to see if there are any segments sitting completely or partially
    in the intertrial space, and trims or removes the segments.
    Arguments:
        df: Pandas dataframe, contains event information. Produced by running
            scripts in preprocessing.
    """
    SEGS_START = ['C1', 'O1']
    SEGS_STOP = ['C2', 'O2']
    TRIALS_START = ['11', '01']
    TRIALS_STOP = ['12', '02']

    trials = df[[code in TRIALS_START + TRIALS_STOP for code in df.Type]]

    # CASE 1: Segment sits wholly in intertrial space. Remove the whole thing.
    bad_idx = []
    for i in range(df.shape[0]):
        if df.iloc[i].Type in SEGS_START and (in_intertrial(df.iloc[i].Latency, trials) and
                                              in_intertrial(df.iloc[i+1].Latency, trials)):
            bad_idx.append(i)
            bad_idx.append(i+1)
    df = df.drop(df.index[bad_idx])
    df = df.set_index(np.arange(0, df.shape[0], 1))

    for i in range(df.shape[0]):
        # CASE 2: Segment head sits inside of the intertrial space. Trim the front of the segment.
        if df.iloc[i].Type in SEGS_START and in_intertrial(df.iloc[i].Latency, trials):
            df.at[i, 'Latency'] = get_next_trial(df.iloc[i].Latency, trials).at[0, 'Latency']

        # CASE 3: Segment tail sits inside of the intertrial space. Trim the end of the segment.
        if df.iloc[i].Type in SEGS_STOP and in_intertrial(df.iloc[i].Latency, trials):
            df.at[i, 'Latency'] = get_previous_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
    return df
