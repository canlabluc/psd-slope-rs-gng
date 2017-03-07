"""
Imports EMSE-exported .evt files and exports them as
a pandas dataframe.
"""

import os
import sys
import glob

import xmltodict
import pandas as pd


def get_filelist(import_path, extension):
    """
    Returns list of file paths from import_path with specified extension.
    """
    filelist = []
    for root, dirs, files in os.walk(import_path):
        filelist += glob.glob(os.path.join(root, '*.' + extension))
        return filelist


def get_event_file(filepath, include_clean_segs=True):
    events = {}
    xmlevents = xmltodict.parse(open(filepath).read())
    xmlevents = xmlevents['EMSE_Event_List']['Event']
    n = 0
    for i in range(len(xmlevents)):
        if include_clean_segs is True and xmlevents[i]['Name'] in ['C', 'O']:
            # First, add front of clean segment
            events[n] = {'type': '', 'latency': 0, 'urevent': 0}
            events[n]['type'] = xmlevents[i]['Name'] + '1'
            events[n]['latency'] = int(xmlevents[i]['Start'])
            events[n]['urevent'] = n

            # Then add the end of the clean segment
            n += 1
            events[n] = {'type': '', 'latency': 0, 'urevent': 0}
            events[n]['type'] = xmlevents[i]['Name'] + '2'
            events[n]['latency'] = int(xmlevents[i]['Stop'])
            events[n]['urevent'] = n
            n += 1
        if len(xmlevents[i]['Name']) == 3 and xmlevents[i]['Name'] not in ['255', '222', '223', '252']:
            events[n] = {'type': '', 'latency': 0, 'urevent': 0}
            events[n]['type'] = xmlevents[i]['Name']
            events[n]['latency'] = int(xmlevents[i]['Start'])
            events[n]['urevent'] = n
            n += 1
    return events


def print_events(files, include_clean_segs=False):
    for file in files:
        evts = get_event_file(file, include_clean_segs)
        print(file.split('/')[-1])
        for i in range(len(evts)):
            print(evts[i]['type'])


def print_evt_information(filepath, print_in_secs=True):
    """
    Prints general information about a specified .evt file.
    Inputs:
        filepath: String, specifies full path to a .evt file.
        print_clean_segs: Boolean, specifies whether to print information regarding
                          the segments that have been marked as being clean.
        print_in_secs: Boolean, specifies whether to print .evt information in seconds,
                       as opposed to timepoints. Assumes 512 sampling rate.
    """
    trials = get_event_file(filepath, include_clean_segs=True)
    # clean_segs = get_event_file(filepath, clean_segs_only=True)
    eyesc_trials = 0
    eyeso_trials = 0
    eyesc_trial_length = 0
    eyeso_trial_length = 0
    for i in range(len(trials)-1):
        # Eyes-closed trial
        if trials[i]['type'][0:2] == '10' and trials[i+1]['type'][0:2] == '20':
            eyesc_trials += 1
            eyesc_trial_length += trials[i+1]['latency'] - trials[i]['latency']
        # Eyes-open trial
        elif trials[i]['type'][0:2] == '11' and trials[i+1]['type'][0:2] == '21':
            eyeso_trials += 1
            eyeso_trial_length += trials[i+1]['latency'] - trials[i]['latency']
    eyesc_trial_length /= eyesc_trials
    eyeso_trial_length /= eyeso_trials
    if print_in_secs:
        eyesc_trial_length /= 512
        eyeso_trial_length /= 512
    print('{0} | eyesc_trials: {1} | eyeso_trials: {2} | eyesc_trial_length: {3} | eyeso_trial_length: {4}'.format(
        filepath.split('/')[-1], eyesc_trials, eyeso_trials, eyesc_trial_length, eyeso_trial_length))


importpath = sys.argv[1]
exportpath = sys.argv[2]

evt_files = get_filelist(importpath, 'evt')
for i in range(len(evt_files)):
    events_dict = get_event_file(evt_files[i], include_clean_segs=True)
    types = [events_dict[j]['type'] for j in range(len(events_dict))]
    latencies = [events_dict[j]['latency'] for j in range(len(events_dict))]

    df = pd.DataFrame(data={'Type': types, 'Latency': latencies})
    df = df[['Type', 'Latency']]
    df.to_csv(exportpath + evt_files[i].split('/')[-1], index=False, sep='\t')
