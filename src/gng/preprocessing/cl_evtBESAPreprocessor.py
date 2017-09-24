#!/usr/local/bin/python3
""" cl_evtBESAPreprocessor
Takes raw BESA-exported evt files and exports them as a cleaner Pandas
dataframe. BESA exports user-defined triggers and event codes in
different columns. This program consolidates them into a single column.

Usage:
    $ ./cl_BESAevtparser.py importpath exportpath

Inputs:
    importpath: Path to directory containing raw evt files.
    exportpath: Path to directory into which to export processed evt files.

Notes:
cl_evtBESAPreprocessor takes files that look like:

    Tmu             Code    TriNo   Comnt
    0               41      2012-10-26T11:27:56.000
    33277           21      0
    2162500         11      0       Pattern 1
    3746484         11      0       Pattern 1
    4793359         11      0       Pattern 1
    ...

And produces files that look like:

    
"""

import os
import sys
import glob
import numpy as np
import pandas as pd

import_path = sys.argv[1]
export_path = sys.argv[2]
clean_file  = True

evt_files = []
for root, dirs, files in os.walk(import_path):
    evt_files += glob.glob(os.path.join(root, '*.evt'))

for i in range(len(evt_files)):
    print('Processing: {}...'.format(evt_files[i].split('/')[-1]), end='')
    df = pd.read_csv(evt_files[i], sep='\t')
    df.columns = ['Tmu', 'Code', 'Trigger', 'Comment']
    df.Trigger = list(map(lambda x: int(x) if len(x) < 3 and x != '-' else -1, list(df.Trigger)))
    df['Latency'] = round((df.Tmu / 10**6) * 512)

    events = []
    latencies = []
    for j in range(df.shape[0]):

        # Grab current event and event latency
        curr_event = df.iloc[j, :]
        latency = curr_event.Latency

        # And match the code to an event
        # CHECK BESA CODES
        if curr_event.Code == 41:
            event = 'INVALID'
        elif curr_event.Code == 11:
            event = 'BLINK1'
        elif curr_event.Code == 21:
            event = 'ARTFCT1'
        elif curr_event.Code == 22:
            event = 'ARTFCT2'
        else:
            # CHECK TRIGGERS
            if curr_event.Trigger == 11:
                event = 'FIXATION'
            elif curr_event.Trigger == 21:
                event = 'GO_PROMPT'
            elif curr_event.Trigger == 22:
                event = 'NOGO_PROMPT'
            elif curr_event.Trigger == 255:
                event = 'START'
            elif curr_event.Trigger == 1:
                event = 'RESPONSE' # Passed when subject responds or at end of 900s if no response
            elif curr_event.Trigger == 32:
                event = 'INCORRECT_RESPONSE'
            elif curr_event.Trigger == 33:
                event = 'NO_RESPONSE'
            else:
                #print('NOT RECOGNIZED {}: \n\tLatency: {}\n\tCode: {}\n\tTrigger: {}\n\tComment: {}'.format(\
                #      evt_files[i].split('/')[-1], event.Latency, event.Code, event.Trigger, event.Comment))
                event = curr_event[1]
        # If the current event is a valid code, add it to the list.
        if event != 'INVALID':
            events.append(event)
            latencies.append(latency)
            if event == 'BLINK1': # Add tail end of blinks
                events.append('BLINK2')
                latencies.append(latency + 205)
    # Construct Pandas dataframe using extracted events and latencies
    df = pd.DataFrame(data={'Event': events, 'Latency': latencies})

    # If the user has specified to return a cleaned file, we remove
    # any trials which contain artifact markers, and remove trials
    # with incorrect or no responses.
    if clean_file:

        ARTIFACT = ['ARTFCT1', 'ARTFCT2', 'BLINK1', 'BLINK2']

        # Removing incorrect trials
        # Find every single INCORRECT_RESPONSE and NO_RESPONSE marker.
        # Then, search backward to find its corresponding FIXATION
        # marker and mark the whole trial as artifact.
        bad_idx = []
        for j in range(df.shape[0]):
            if df.iloc[j].Event in ['INCORRECT_RESPONSE', 'NO_RESPONSE']:
                prior = df[df.Latency < df.iloc[j].Latency]
                fixation_idx = prior[prior.Event.isin(['FIXATION'])].tail(1).index[0]
                response_idx = prior[prior.Event.isin(['RESPONSE'])].tail(1).index[0]

                # Replace the FIXATION and INCORR_RESP/NO_RESP markers
                # with artifact markers, and throw away the RESPONSE
                # marker between them.
                df.set_value(fixation_idx, 'Event', 'ARTFCT1')
                df.set_value(j,            'Event', 'ARTFCT2')
                # df.at[fixation_idx, 'Event'] = 'ARTFCT1'
                # df.at[j,            'Event'] = 'ARTFCT2'
                bad_idx.append(response_idx)

                # print('fixation_idx {}'.format(fixation_idx))
                # print('response_idx {}'.format(response_idx))
                # print('incorr/no    {}'.format(j))
        df = df.drop(df.index[bad_idx])
        df = df.set_index(np.arange(0, df.shape[0], 1))

        # Removing correct-response trials that contain artifact
        # For each FIXATION marker, check to see if that trial contains
        # artifact markers. If it does, we throw it out as above.
        for j in range(df.shape[0]):
            if df.iloc[j].Event == 'FIXATION':
                # Grab the corresponding RESPONSE marker's row index in
                # the dataframe.
                upcoming = df[df.Latency > df.iloc[j].Latency]
                response_idx = upcoming[upcoming.Event.isin(['RESPONSE'])].head(1).index[0]

                # Check if there are any ARTFCT markers between FIXATION
                # and RESPONSE
                markers = [df.iloc[k].Event for k in range(j+1, response_idx)]
                if any(m in markers for m in ARTIFACT):
                    # If there are, replace FIXATION and RESPONSE
                    # markers with artifact markers.
                    df.at[j,            'Event'] = 'ARTFCT1'
                    df.at[response_idx, 'Event'] = 'ARTFCT2'

    subj_name = evt_files[i].split('/')[-1][:-4]
    print(' Done.')
    df.to_csv(export_path + '/' + subj_name + '.evt', sep='\t', index=False)
