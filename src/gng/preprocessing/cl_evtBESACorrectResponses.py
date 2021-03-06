#!/usr/local/bin/python3
""" cl_evtBESACorrectResponses
Constructs evt file that contains only artifact-free correct responses in the
Go-NoGo task. We obtain the Go-NoGo triggers from preprocessing the evt files
using cl_evtBESAPreprocessor, and the clean segments from cl_evtBESACleanSegments.

Usage:
    $ ./cl_evtBESACorrectResponses.py importpath_preprocessed importpath_clean exportpath

Inputs:
    importpath_preprocessed: Directory containing evt files that have been preprocessed
                             using cl_evtBESAPreprocessor.
    importpath_clean: Directory containing evt files that have been preprocessed using
                      cl_evtBESAPreprocessor and then through cl_evtBESACleanSegments.
    exportpath: Directory in which to export processed evt files.

Notes:
    cl_evtBESACorrectResponses takes file that look like:

    Latency Trigger
    17.0    ARTFCT1
    1107.0  BLINK1
    1918.0  BLINK1
    2454.0  BLINK1
    ...

    Latency Trigger
    207252.0        C1
    210038.0        C2
    210243.0        C1
    ...

    and produces a file that looks like:

    Latency Trigger
    209897.0        NOGO_PROMPT
    210426.0        RESPONSE
    214680.0        NOGO_PROMPT
"""

import os
import sys
import glob

import numpy as np
import pandas as pd

importpath_conso = sys.argv[1]
importpath_clean = sys.argv[2]
exportpath       = sys.argv[3]

evt_files = []
for root, dirs, files in os.walk(importpath_conso):
    evt_files += glob.glob(os.path.join(root, '*.evt'))

for i in range(len(evt_files)):
    # Open the file, retaining only prompt and response-related data.
    subj_name = evt_files[i].split('/')[-1][:-4]
    df      = pd.read_csv(importpath_conso + subj_name + '.evt', sep='\t')
    cleandf = pd.read_csv(importpath_clean + subj_name + '.evt', sep='\t')
    df = df[df.Trigger.isin(['GO_PROMPT', 'NOGO_PROMPT', 'RESPONSE', \
                             'INCORRECT_RESPONSE', 'NO_RESPONSE'])]
    df = df.set_index(np.arange(0, df.shape[0], 1))

    # Filter out incorrect responses.
    bad_idx = []
    for i in range(df.shape[0] - 3):
        if df.iloc[i].Trigger == 'GO_PROMPT' and (df.iloc[i+2].Trigger == 'NO_RESPONSE' or
                                                  df.iloc[i+2].Trigger == 'INCORRECT_RESPONSE'):
            bad_idx.append(i); bad_idx.append(i+1); bad_idx.append(i+2)
        if df.iloc[i].Trigger == 'NOGO_PROMPT' and df.iloc[i+2].Trigger == 'INCORRECT_RESPONSE':
            bad_idx.append(i); bad_idx.append(i+1); bad_idx.append(i+2)
    df = df.drop(df.index[bad_idx])
    df = df.set_index(np.arange(0, df.shape[0], 1))

    # Filter out prompts/responses which are not contained in clean segments.
    bad_idx = []
    for i in range(df.shape[0] - 3):
        if df.iloc[i].Trigger == 'GO_PROMPT' or df.iloc[i].Trigger == 'NOGO_PROMPT':
            # 1. Check that it's inside of a clean segment
            latency = df.iloc[i].Latency
            if not cleandf[cleandf.Latency < latency].tail(1).Trigger.all() == 'C1':
                bad_idx.append(i); bad_idx.append(i+1)
    df = df.drop(df.index[bad_idx])
    df = df.set_index(np.arange(0, df.shape[0], 1))

    # Filter out prompts/responses which do not have 600ms of post-stimulus clean space.
    bad_idx = []
    for i in range(df.shape[0] - 3):
        if df.iloc[i].Trigger == 'GO_PROMPT' or df.iloc[i].Trigger == 'NOGO_PROMPT':
            latency = df.iloc[i].Latency
            end_clean_seg = cleandf[cleandf.Latency > latency].iloc[0].Latency
            if not end_clean_seg - latency < 308: # 0.6s * 512 = 307.2 points
                bad_idx.append(i); bad_idx.append(i+1)
    df = df.drop(df.index[bad_idx])
    df = df.set_index(np.arange(0, df.shape[0], 1))

    df.to_csv(exportpath + '/' + subj_name + '.evt', sep='\t', index=False)
