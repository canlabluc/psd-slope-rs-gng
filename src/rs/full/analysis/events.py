import numpy as np
import pandas as pd

"""
Set of functions for cleaning up event dataframes.
"""

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
