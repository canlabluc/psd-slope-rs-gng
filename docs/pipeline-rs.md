# preprocessing and data analysis pipeline for resting state data

This file specifies the pipeline utilized to process resting-state data for studies 127 and 118. Generally, processing occurs in the following way:

First, preprocess .evt files in order to clean them up for later use. If this has already been done, then skip this step.
```bash
$ python rm_irrelevant_xml.py
$ python transform_data_to_df.py
$ python rm_entire_trial_segs.py
$ python rm_intertrial_segs.py
``` 

From here we pre-preprocess files in MATLAB, utilizing the `cl_preprocessingoriginal.m` script. This calls other functions in the CAN Lab toolbox in order to import raw data, import event markers, re-reference data, apply a band-pass filter, construct channel clusters if we're processing sensor-level data, and lastly export data to .mat format for importing into Python. 

