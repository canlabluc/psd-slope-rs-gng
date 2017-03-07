# resting state protocol for studies 127 and 118

This file specifies the protocol under which the resting state data was taken for studies 127 and 118. For both datasets, partipants had EEG recorded during two conditions: eyes closed and eyes open. The pattern and duration of the eyes-closed and eyes-open trials differ slightly between the two datasets. Event-related information, exported from EMSE, is found in the .evt files.

# 127 - Older adults
Each trial has a duration of 30 seconds.
Participants undergo a total of 8 trials. 4 are eyes-closed, 4 are eyes-open. Thus, there's 2 minutes of eyes-closed data and 2 minutes of eyes-open data.
Time between trials is variable, typically between 5 - 15 seconds.

Trials are organized into one of the two following sequences, where `C` is eyes-closed and `O` is eyes-open:
```
COOC OCCO
```
```
OCCO COOC
```

# 118 - Younger adults
Each trial has a duration of 60 seconds.
Participants undergo a total of 8 trials. 4 are eyes-closed, 4 are eyes-open.
Thus, there's 4 minutes of eyes-closed data and 4 minutes of eyes-open data.
Time between trials is variable, typically between 5 - 20 seconds.

Trials are organized into the same sequences as the older adults:
```
COOC OCCO
```
```
OCCO COOC
```

# How event-related files are processed
Files specifying events in the EEG recording are treated in the same general manner for both resting-state and task-related data, in that the files are exported from either EMSE (resting-state) or BESA (task-related), cleaned up either manually or through scripts, and then fed into the pipeline.

Resting-state data is exported from EMSE in an xml format. To check data integrity, a set of scripts were written in Python. These scripts do the following:

1. Use the `xmltodict` Python library to export the xml files to Python dictionaries. From here, the dictionary is translated into a Pandas dataframe with the following format:
```
	Latency Type
0 0				[seg]
1 1244		222
2 2484		101
...
```
2. The Pandas dataframe makes it a bit easier to work with the event files. We then utilize this file to check a few things:
- Check that files have the correct number of `C` and `O` trials, of the right length.
- Check that segments that have been marked clean in the .evt file are sitting inside of a trial. In other words, that there are no segments of data marked as clean inside intertrial periods.

# Issues with data
The following subjects had strange .evt files that required further inspection.

### 112118266.evt
Contains only 4 full trials: two closed and two open, in the `OCCO` sequence. It seems the event file got cut off in the middle of the 5th trial.

### 120127123.evt
Contains only 7 full trials. This subject is missing the very first `O` trial.

### 120127148.evt
Contains only 7 full trials. This subject is missing the beginning marker for the first trial, an `O` trial. This was fixed in the `evt-data` repository. The beginning marker was added manually at time-point 0.

### 120127170.evt
Contains 8 full trials, but the first trial (an eyes-closed trial) contains no marked clean segments.
