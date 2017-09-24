"""
This script cleans up EMSE-exported .evt files, through the use of
a few different functions.

rm_irrelevant_xml removes unneeded markers from the xml files.
These are the following: '[seg]', '222', '252', '223', '255', as well as
the unneeded headings at the top of the files, such as <bUseQID>.

Once these are clean, transform_xml_to_df imports the files and produces
tab-separated csv files containing the event information, using Pandas.
"""

import os
import sys
import glob
import getopt

import xmltodict
import numpy as np
import pandas as pd

################################################################################

# Auxilliary functions

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


def print_seg_code_information(df, fname, i, trials, error_type):
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


def get_cmdline_params(params):
    help_msg = '\tevt_preprocessing.py -i <import_dir> -o <export_dir>\n\n'
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
    except getopt.GetoptError:
        print('Error: Bad input. To run:\n')
        print(help_msg)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Help:\n')
            print(help_msg)
            sys.exit(2)
        elif opt == '-i':
            params['import_path'] = arg
        elif opt == '-o':
            params['export_path'] = arg
    return params


################################################################################

# Preprocessing functions

def rm_irrelevant_xml(import_path, export_path):
    """
    Removes unneeded xml tags from EMSE .evt files.

    Parameters
    ----------
    import_path : str
        Absolute path to EMSE-exported .evt files.

    export_path : str
        Absolute path to directory in which to save cleaned up EMSE
        .evt files.
    """
    files = get_filelist(import_path, 'evt')
    for file in files:
        lines = open(file, 'r').readlines()
        lines[1] = ''
        lines[2] = ''
        lines[3] = ''
        for i in range(len(lines)):
            if ('<Name>[seg]</Name>' in lines[i] or
               '<Name>222</Name>' in lines[i] or
               '<Name>223</Name>' in lines[i] or
               '<Name>252</Name>' in lines[i] or
               '<Name>255</Name>' in lines[i]):
               lines[i] = ''
        new = open(export_path + file.split('/')[-1], 'w')
        new.writelines(lines)


def transform_xml_to_df(import_path, export_path):
    """
    Transforms EMSE .evt files from xml to tab-separated Pandas
    dataframes. Note that input to these should be cleaned up
    with remove_irrelevant_xml.

    Parameters
    ----------
    import_path : str
        Absolute path to EMSE-exported .evt files in xml format,
        cleaned through remove_irrelevant_xml.

    export_path : str
        Absolute path to directory in which to save tab-separated
        Pandas dataframe .evt files.
    """

    evt_files = get_filelist(import_path, 'evt')
    for i in range(len(evt_files)):
        events_dict = get_event_file(evt_files[i], include_clean_segs=True)
        types = [events_dict[j]['type'] for j in range(len(events_dict))]
        latencies = [events_dict[j]['latency'] for j in range(len(events_dict))]

        df = pd.DataFrame(data={'Type': types, 'Latency': latencies})
        df = df[['Type', 'Latency']]
        df.to_csv(export_path + evt_files[i].split('/')[-1], index=False, sep='\t')


def rm_entire_trial_segs(import_path, export_path):
    TYPES_START = ['11', '10']
    TYPES_STOP  = ['21', '20']
    SEGS_START  = ['C1', 'O1']
    SEGS_STOP   = ['C2', 'O2']
    files = get_filelist(import_path, 'evt')
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

        df.to_csv(export_path + f.split('/')[-1], sep='\t', index=False)


def rm_intertrial_segs(import_path, export_path):
    SEGS_START = ['C1', 'O1']
    SEGS_STOP = ['C2', 'O2']
    TRIALS_START = ['11', '10']
    TRIALS_STOP = ['21', '20']
    files = get_filelist(import_path, 'evt')
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
                # print_seg_code_information(df, i, trials, 'INTERTRIAL')
                bad_idx.append(i)
                bad_idx.append(i+1)
                found_bad_segs = True
        df = df.drop(df.index[bad_idx])
        df = df.set_index(np.arange(0, df.shape[0], 1))

        # Trim clean segments that sit partially in intertrial space.
        for i in range(df.shape[0]):

            # CASE 2: Segment head sits inside of the intertrial space. Trim the front of the segment.
            if df.iloc[i].Type[0:2] in SEGS_START and in_intertrial(df.iloc[i].Latency, trials):
                # print_seg_code_information(df, fname, i, trials, 'TRIM HEAD')
                df.at[i, 'Latency'] = get_next_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
                found_bad_segs = True

            # CASE 3: Segment tail sits inside of the intertrial space. Trim the end of the segment.
            if df.iloc[i].Type[0:2] in SEGS_STOP and in_intertrial(df.iloc[i].Latency, trials):
                # print_seg_code_information(df, fname, i, trials, 'TRIM TAIL')
                df.at[i, 'Latency'] = get_previous_trial(df.iloc[i].Latency, trials).at[0, 'Latency']
                found_bad_segs = True

        if found_bad_segs:
            print()

        df.to_csv(export_path + fname, sep='\t', index=False)


## Parameters ##################################################################


def main(argv):
    """
    Preprocessing .evt files.

    Parameters
    ----------
    import_path : str
        Absolute path to EMSE-exported .evt files in xml format,
        cleaned through remove_irrelevant_xml.

    export_path : str
        Absolute path to directory in which to save tab-separated
        Pandas dataframe, preprocessed .evt files.
    """

    params = dict()
    params = get_cmdline_params(params)

    print('Processing .evt files... ', end='')
    rm_irrelevant_xml(params['import_path'], params['export_path'])
    transform_xml_to_df(params['export_path'], params['export_path'])
    rm_entire_trial_segs(params['export_path'], params['export_path'])
    # rm_intertrial_segs(params['export_path'], params['export_path']) # TODO: Fix

    # Subject 112118266 contains an extra erroneous trial which we
    # simply remove from the processed .evt file:
    files = get_filelist(params['export_path'], 'evt')
    for i in range(len(files)):
        if files[i].split('/')[-1] == '112118266.evt':
            with open(files[i]) as bad_file:
                lines = bad_file.readlines()
            w = open(files[i], 'w')
            w.writelines([event for event in lines[:-1]])
            w.close()

    print('Done.')


if __name__ == '__main__':
    main(sys.argv)