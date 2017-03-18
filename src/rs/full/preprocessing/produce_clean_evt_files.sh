#!/bin/bash

raw_evt_files='../../../../data/rs/full/evt/raw/'
clean_evt_files='../../../../data/rs/full/evt/clean/'

python rm_irrelevant_xml.py $raw_evt_files $clean_evt_files
python transform_data_to_df.py $clean_evt_files $clean_evt_files
python rm_entire_trial_segs.py $clean_evt_files $clean_evt_files
python rm_intertrial_segs.py $clean_evt_files $clean_evt_files
