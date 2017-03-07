# log

## March 6, 2017
#### Restructured src/rs/full, and updated evt data.
src/rs/full has been restructured into preprocessing and analysis folders.

The resting-state evt files we were using for older adults had some problems. See the [evt-data](https://github.com/canlabluc/evt-data) project for more information. This might be a confound for the PSD analysis because some trials in the recording (primarily eyes-closed trials) were being wholly marked as clean data. EMG artifact, which I assume trials sometimes contained, is generally characterized by high-frequency oscillations, which could account for the flatter PSD slope in the older adults.

Additionally, some of the segments that were marked as being clean data were sitting either completely in the intertrial period or partially inside of it. The preprocessing scripts fix this.

In order to address these issues, I'm running all of our full recording data through the analysis again, but this time with the cleaned up evt files. We'll be able to see how much this impacts:

1. Number of windows able to be extracted for PSD computation (there will likely be considerably fewer extractable windows)
2. The difference in neural noise between older and younger adults

Here are the steps to do so:

Producing the clean evt files:
```
$ bash src/rs/full/preprocessing/produce_clean_evt_files.sh
```

Then, we change the event information for all of the source models and the sensor-level data, exporting them to new folders with the `MagCleanEvtFiltCAR` designation:

```matlab
% Default Mode Network
mkdir 'data/rs/full/source-dmn/MagCleanEvtFiltCAR-set'
mkdir 'data/rs/full/source-dmn/MagCleanEvtFiltCAR-mat'
cl_modifyevents('data/rs/full/source-dmn/MagEvtFiltCAR-set/',...
             'data/rs/full/evt/clean/',...
             'data/rs/full/source-dmn/',...
             {'C', 'O'});
set_mat_converter('data/rs/full/source-dmn/MagCleanEvtFiltCAR-set/',...
                  'data/rs/full/source-dmn/MagEvtFiltCAR-mat/');

% Frontal
mkdir 'data/rs/full/source-frontal/MagCleanEvtFiltCAR-set'
mkdir 'data/rs/full/source-frontal/MagCleanEvtFiltCAR-mat'
cl_modifyevents('data/rs/full/source-frontal/MagEvtFiltCAR-set/',...
             'data/rs/full/evt/clean/',...
             'data/rs/full/source-frontal/',...
             {'C', 'O'});
set_mat_converter('data/rs/full/source-frontal/MagCleanEvtFiltCAR-set/',...
                  'data/rs/full/source-frontal/MagCleanEvtFiltCAR-mat/');

% Dorsal
mkdir 'data/rs/full/source-dorsal/MagCleanEvtFiltCAR-set'
mkdir 'data/rs/full/source-dorsal/MagCleanEvtFiltCAR-mat'
cl_modifyevents('data/rs/full/source-dorsal/MagEvtFiltCAR-set/',...
             'data/rs/full/evt/clean/',...
             'data/rs/full/source-dorsal/',...
             {'C', 'O'});
set_mat_converter('data/rs/full/source-dorsal/MagCleanEvtFiltCAR-set/',...
                  'data/rs/full/source-dorsal/MagCleanEvtFiltCAR-mat/');

% Ventral
mkdir 'data/rs/full/source-ventral/MagCleanEvtFiltCAR-set'
mkdir 'data/rs/full/source-ventral/MagCleanEvtFiltCAR-mat'
cl_modifyevents('data/rs/full/source-ventral/MagEvtFiltCAR-set/',...
             'data/rs/full/evt/clean/',...
             'data/rs/full/source-ventral/',...
             {'C', 'O'});
set_mat_converter('data/rs/full/source-ventral/MagCleanEvtFiltCAR-set/',...
                  'data/rs/full/source-ventral/MagCleanEvtFiltCAR-mat/');

% Sensor-level data
mkdir 'data/rs/full/original/ExclFiltCleanEvtCARClust-set'
mkdir 'data/rs/full/original/ExclFiltCleanEvtCARClust-mat'
cl_modifyevents('data/rs/full/original/ExclFiltCleanEvtCARClust-set',...
             'data/rs/full/evt/clean/',...
             'data/rs/full/original/ExclFiltCleanEvtCARClust-set/',...
             {'C', 'O'});
set_mat_converter('data/rs/full/original/ExclFiltCleanEvtCARClust-set/',...
                  'data/rs/full/original/ExclFiltCleanEvtCARClust-mat/');
```

Then we run `spectral_slopes.py` on the new data:
```bash
$ python src/rs/full/analysis/spectral_slopes.py -m source-level -i data/rs/full/original/ExclFiltCleanEvtCARClust-mat/ -o data/runs/
$ python src/rs/full/analysis/spectral_slopes.py -m dmn -i data/rs/full/source-dmn/MagCleanEvtFiltCAR-mat/ -o data/runs/
$ python src/rs/full/analysis/spectral_slopes.py -m frontal -i data/rs/full/source-frontal/MagCleanEvtFiltCAR-mat/ -o data/runs/
$ python src/rs/full/analysis/spectral_slopes.py -m ventral -i data/rs/full/source-ventral/MagCleanEvtFiltCAR-mat/ -o data/runs/
$ python src/rs/full/analysis/spectral_slopes.py -m dorsal -i data/rs/full/source-dorsal/MagCleanEvtFiltCAR-mat/ -o data/runs/
```
