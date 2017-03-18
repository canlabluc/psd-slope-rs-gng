"""
This file removes segments which indicate when a trial begins
and ends, but have erroneously been marked as being clean segment
indicators.

This operates on df-tranformed evt files. Run transform_data_to_df.py
prior to running this script.
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


importpath = sys.argv[1]
exportpath = sys.argv[2]

TYPES_START = ['11', '10']
TYPES_STOP = ['21', '20']
SEGS_START = ['C1', 'O1']
SEGS_STOP = ['C2', 'O2']
files = get_filelist(importpath, 'evt')
for f in files:

    df = pd.read_csv(f, sep='\t')

    bad_idx = []
    for i in range(df.shape[0]):
        # Check whether we're currently in intertrial or intratrial space.
        if df.iloc[i].Type[0:2] in TYPES_START:
            intertrial = False
        elif df.iloc[i].Type[0:2] in TYPES_STOP:
            intertrial = True

        # If we find a segment marked as clean in the intertrial period, mark
        # it as bad and remove it from the .evt file.
        if (df.iloc[i].Type[0:2] in SEGS_START and
           intertrial is True and
           (df.iloc[i+1].Latency - df.iloc[i].Latency)/512 > 29):
            bad_idx.append(i)
            bad_idx.append(i+1)
    df = df.drop(bad_idx)
    df = df.set_index(np.arange(0, df.shape[0], 1))

    df.to_csv(exportpath + f.split('/')[-1], sep='\t', index=False)
