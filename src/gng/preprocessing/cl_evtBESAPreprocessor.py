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
import pandas as pd

import_path = sys.argv[1]
export_path = sys.argv[2]

print(import_path)
print(export_path)

evt_files = []
for root, dirs, files in os.walk(import_path):
    evt_files += glob.glob(os.path.join(root, '*.evt'))

for i in range(len(evt_files)):
    df = pd.read_csv(evt_files[i], sep='\t')
    df.columns = ['Tmu', 'Code', 'Trigger', 'Comment']
    df.Trigger = list(map(lambda x: int(x) if len(x) < 3 and x != '-' else -1, list(df.Trigger)))
    df['Latency'] = round((df.Tmu / 10**6) * 512)

    events = []
    latencies = []
    for j in range(df.shape[0]):

        # Grab current event and event latency
        curr_event = df.iloc[j, :]
        latency = event.Latency

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
                event = event[1]
        if event != 'INVALID':
            events.append(event)
            latencies.append(latency)
            if event == 'BLINK1': # Add tail end of blinks
                events.append('BLINK2')
                latencies.append(latency + 205)

    df_export = pd.DataFrame(data={'Event': events, 'Latency': latencies})
    subj_name = evt_files[i].split('/')[-1][:-4]
    print("Processed: " + export_path + '/' + subj_name + '.evt')
    df_export.to_csv(export_path + '/' + subj_name + '.evt', sep='\t', index=False)
