# data pipeline for go-nogo data

This file specifies the pipeline utilized to process gng data for the older and younger adults. Generally, processing occurs in the following way:

#### preprocessing
Most of the preprocessing occurs in MATLAB, using the CAN Lab EEGLAB plugin. This is handled by `cl_preprocessinggng.m`. Simply modify the paths prior to running it.

After running `cl_preprocessinggng.m`, we clean up .evt files in order to use them in the analysis (using the Python scripts in `src/gng/preprocessing/`):
```bash
$ bash produce_clean_evt_files.sh
```

#### spectral slopes analysis
With EEG data preprocessed and .evt files ready for use, 




