#!/bin/bash

# For running check_YA_trials.py

# Define analysis, importpaths, exportpaths
analysis='src/rs/full/analysis/check_YA_trials.py'
import_dir_sensor='data/rs/full/original/ExclFiltCleanEvtCARClust-mat'
import_dir_dmn='data/rs/full/source-dmn/MagCleanEvtFiltCAR-mat'
import_dir_frontal='data/rs/full/source-frontal/MagCleanEvtFiltCAR-mat'
import_dir_dorsal='data/rs/full/source-dorsal/MagCleanEvtFiltCAR-mat'
import_dir_ventral='data/rs/full/source-ventral/MagCleanEvtFiltCAR-mat'
export_dir='data/runs/2017-03-13/'

# Take PSDs of first half of trials.
python $analysis -m dmn     -i $import_dir_dmn     -o $export_dir -c match_OA
python $analysis -m frontal -i $import_dir_frontal -o $export_dir -c match_OA
python $analysis -m dorsal  -i $import_dir_dorsal  -o $export_dir -c match_OA
python $analysis -m ventral -i $import_dir_ventral -o $export_dir -c match_OA
python $analysis -m sensor-level -i $import_dir_sensor -o $export_dir -c match_OA

# Second half of trials.
python $analysis -m dmn     -i $import_dir_dmn     -o $export_dir -c latter_half
python $analysis -m frontal -i $import_dir_frontal -o $export_dir -c latter_half
python $analysis -m dorsal  -i $import_dir_dorsal  -o $export_dir -c latter_half
python $analysis -m ventral -i $import_dir_ventral -o $export_dir -c latter_half
python $analysis -m sensor-level -i $import_dir_sensor -o $export_dir -c latter_half

