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

# General pipeline for computing PSD slopes
TODO

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

### 120127122.evt
Contained three events called "Closed2" which seem to span lengths of time much larger than typical clean segments. It seems to likely be marking the trial. Not sure why these are here, but they've been removed.

### 120127123.evt
Contains only 7 full trials. This subject is missing the very first `O` trial.

### 120127148.evt
Contains only 7 full trials. This subject is missing the beginning marker for the first trial, an `O` trial.

### 120127170.evt
Contains 8 full trials. There is an extra eyes-closed trial in the beginning that was input erroneously. There are blinks during the trial, indicating that the subject had their eyes open. This trial was removed from .evt file.

### Data integrity: segments that were found in intertrial period (bad)
```
In [8]: for file in files:
   ...:     name = file.split('/')[-1]
   ...:     print(name)
   ...:     evts = get_event_file(file, True)
   ...:     intertrial = True
   ...:     for i in range(len(evts)):
   ...:         if evts[i]['type'][0:2] in ['21', '20']:
   ...:             intertrial = True
   ...:         elif evts[i]['type'][0:2] in ['11', '10']:
   ...:             intertrial = False
   ...:         
   ...:         if intertrial == True and evts[i]['type'] in ['C1', 'O1']:
   ...:             print('\t {} | Start: {} | Stop: {}'.format(evts[i]['type'], evts[i]['latency'], evts[i+1]['latency']))
   ...:             
   ...:         
1121181181.evt
1121181183.evt
1121181218.evt
1121181262.evt
1121181286.evt
112118131.evt
1121181334.evt
112118135.evt
1121181393.evt
1121181418.evt
1121181424.evt
1121181428.evt
1121181510.evt
1121181517.evt
1121181575.evt
112118167.evt
112118204.evt

112118257.evt
112118266.evt
112118334.evt
112118373.evt
112118416.evt
112118463.evt
112118468.evt
112118475.evt
112118479.evt
112118521.evt
112118526.evt
112118576.evt
112118578.evt
112118587.evt
112118642.evt
112118723.evt
112118761.evt
112118762.evt
112118785.evt
120127101.evt
120127102.evt
	 C1 | Start: 100898 | Stop: 103253
120127103.evt
120127104.evt
	 O1 | Start: 23750 | Stop: 27229
	 C1 | Start: 126457 | Stop: 141817
120127105.evt
	 O1 | Start: 75025 | Stop: 80014
120127106.evt
120127107.evt
	 O1 | Start: 13926 | Stop: 29696
	 O1 | Start: 82534 | Stop: 98304
120127108.evt
	 C1 | Start: 91884 | Stop: 95366
	 C1 | Start: 113942 | Stop: 125922
	 C1 | Start: 177051 | Stop: 186677
120127109.evt
120127110.evt
	 O1 | Start: 28927 | Stop: 37938
	 O1 | Start: 148877 | Stop: 164237
120127111.evt
120127112.evt
	 O1 | Start: 16179 | Stop: 22835
120127113.evt
	 C1 | Start: 96901 | Stop: 109189
120127114.evt
	 O1 | Start: 165252 | Stop: 168836
120127115.evt
120127116.evt
	 O1 | Start: 30819 | Stop: 36349
120127117.evt
120127118.evt
120127119.evt
	 O1 | Start: 36579 | Stop: 51939
	 C1 | Start: 127857 | Stop: 136151
120127120.evt
	 O1 | Start: 115998 | Stop: 131358
	 O1 | Start: 156170 | Stop: 171530
	 O1 | Start: 175626 | Stop: 176957
120127121.evt
120127122.evt
120127123.evt
	 O1 | Start: 45102 | Stop: 60564
	 O1 | Start: 87732 | Stop: 102989
120127124.evt
120127125.evt
120127128.evt
	 O1 | Start: 30655 | Stop: 36819
	 O1 | Start: 153837 | Stop: 169231
120127130.evt
	 C1 | Start: 70543 | Stop: 85907
	 O1 | Start: 141305 | Stop: 144774
	 C1 | Start: 187465 | Stop: 189893
120127131.evt
	 C1 | Start: 174168 | Stop: 189559
120127132.evt
	 O1 | Start: 13679 | Stop: 15948
	 O1 | Start: 130964 | Stop: 146339
120127133.evt
	 O1 | Start: 29388 | Stop: 32428
120127134.evt
	 O1 | Start: 28925 | Stop: 33649
	 C1 | Start: 112653 | Stop: 128016
120127135.evt
	 C1 | Start: 60690 | Stop: 76062
	 C1 | Start: 87530 | Stop: 102902
120127137.evt
120127138.evt
	 C1 | Start: 68815 | Stop: 84189
	 C1 | Start: 108996 | Stop: 124370
120127139.evt
	 C1 | Start: 38878 | Stop: 54275
	 O1 | Start: 79081 | Stop: 84808
	 C1 | Start: 100102 | Stop: 115477
	 O1 | Start: 119892 | Stop: 124253
120127140.evt
	 O1 | Start: 125769 | Stop: 128495
	 C1 | Start: 165199 | Stop: 180563
120127142.evt
	 O1 | Start: 41685 | Stop: 57050
	 C1 | Start: 63971 | Stop: 79335
	 C1 | Start: 85647 | Stop: 89076
	 O1 | Start: 150627 | Stop: 166013
	 O1 | Start: 171616 | Stop: 186980
	 C1 | Start: 192943 | Stop: 208307
120127144.evt
	 O1 | Start: 26875 | Stop: 29930
	 C1 | Start: 48507 | Stop: 63894
	 O1 | Start: 93513 | Stop: 99426
	 C1 | Start: 118251 | Stop: 133616
	 O1 | Start: 141813 | Stop: 145638
120127145.evt
	 C1 | Start: 101013 | Stop: 116401
	 C1 | Start: 222191 | Stop: 237578
120127146.evt
	 O1 | Start: 114016 | Stop: 115282
120127147.evt
	 O1 | Start: 36532 | Stop: 40849
	 C1 | Start: 61662 | Stop: 74151
	 C1 | Start: 86487 | Stop: 88536
	 O1 | Start: 174327 | Stop: 175473
	 C1 | Start: 195088 | Stop: 207180
120127148.evt
	 O1 | Start: 0 | Stop: 1079
120127149.evt
	 C1 | Start: 39594 | Stop: 49859
120127151.evt
	 C1 | Start: 47219 | Stop: 62549
	 O1 | Start: 90873 | Stop: 103185
120127153.evt
	 C1 | Start: 26651 | Stop: 40857
120127154.evt
120127155.evt
120127156.evt
	 O1 | Start: 40293 | Stop: 42980
	 O1 | Start: 98651 | Stop: 100743
	 C1 | Start: 119478 | Stop: 122716
	 C1 | Start: 177134 | Stop: 187883
120127157.evt
	 C1 | Start: 49446 | Stop: 64821
	 C1 | Start: 115099 | Stop: 130451
	 O1 | Start: 139112 | Stop: 144729
120127158.evt
	 C1 | Start: 90105 | Stop: 93448
	 C1 | Start: 110466 | Stop: 116612
	 O1 | Start: 196686 | Stop: 198411
120127159.evt
	 C1 | Start: 205577 | Stop: 206421
	 C1 | Start: 206484 | Stop: 207426
	 C1 | Start: 208044 | Stop: 208519
	 C1 | Start: 209023 | Stop: 209302
	 C1 | Start: 209604 | Stop: 210101
	 C1 | Start: 210606 | Stop: 211006
120127160.evt
120127161.evt
120127162.evt
120127163.evt
120127164.evt
120127165.evt
120127166.evt
120127167.evt
120127168.evt
120127169.evt
	 O1 | Start: 171461 | Stop: 174260
120127170.evt
	 C1 | Start: 66557 | Stop: 81899
	 C1 | Start: 211389 | Stop: 214687
	 O1 | Start: 231100 | Stop: 238746
```