#raw_evt_files='../../../../../data/gng/evt/raw/'
#df_evt_files='../../../../../data/gng/evt/df-transformed/'
#clean_evt_files='../../../../../data/gng/evt/clean/'

raw_evt_files='/Users/jorge/evts/raw/'
df_evt_files='/Users/jorge/evts/clean/'
clean_evt_files='/Users/jorge/evts/clean/'

python cl_evtBESAPreprocessor.py $raw_evt_files $df_evt_files
