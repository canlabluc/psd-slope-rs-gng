#!/bin/bash
# Handles preprocessing of event-related .evt files

preprocessing="rs"  # Options: rs for resting state, gng for Go-NoGo
raw_evt_files="data/raw-evt/"
clean_evt_files="data/clean-evt/"

if [ "$preprocessing" == "rs" ]; then
  python src/rs/full/preprocessing/rm_irrelevant_xml.py    $raw_evt_files $clean_evt_files
  python src/rs/full/preprocessing/transform_data_to_df.py $clean_evt_files $clean_evt_files
  python src/rs/full/preprocessing/rm_entire_trial_segs.py $clean_evt_files $clean_evt_files
  python src/rs/full/preprocessing/rm_intertrial_segs.py   $clean_evt_files $clean_evt_files
elif [ "$preprocessing" == "gng" ]; then
  python src/gng/preprocessing/cl_evtBESAPreprocessor.py  $raw_evt_files $clean_evt_files
  python src/gng/preprocessing/cl_evtBESACleanSegments.py $raw_evt_files $clean_evt_files
else
  echo "ERROR: Invalid type of preprocessing selected. Use 'gng' or 'rs'."
fi

