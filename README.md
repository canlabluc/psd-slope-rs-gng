# psd-slope

This project organizes the efforts to compute the differences in PSD slope between aMCIs, SAs, and their controls, for both task-related and resting-state data. PSDs are computed at both the sensor and source level (source modeling is done through BESA).

## project layout

```
.
├── data
│   ├── auxilliary
│   ├── gng
│   ├── rs
│   │   ├── 20s
│   │   └── full
│   │       ├── evt
│   │       ├── original
│   │       ├── source-dmn
│   │       ├── source-dorsal
│   │       ├── source-frontal
│   │       └── source-ventral
│   ├── runs
│   └── topographic-corr
│       └── excluded
├── docs
├── figures
├── papers
├── results
├── runs
├── src
│   ├── gng
│   └── rs
│       ├── 20s
│       ├── full-original
│       └── full-src-model
└── work
```