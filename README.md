# psd-slope

This project organizes efforts to compute the differences in PSD slope between aMCIs, SAs, and their controls, for both task-related (Go-NoGo, referred to as `gng`) and resting-state data (referred to as `rs`). PSDs are computed at both the sensor and source level (source modeling is done through BESA).

For instructions on running a specific analysis, see **docs**.

## project layout

The project is organized into a few different folders:
- **data** organizes project data into four directories:
    - **auxilliary:** Contains csv files detailing participant behavioral and non-EEG data.
    - **gng:** Contains EEG data for the Go-NoGo task.
    - **rs:** Contains resting-state EEG data. These are further subdivided into 20-second segments and full recordings. For the full recordings, we have the **original** 66-channel recordings, as well as BESA-exported source models.
    - **runs:** Separate, time-stamped analyses which result from running analyses located in `src` are stored here.
- **docs** contains the steps and protocols for each analysis.
- **figures** contains figures produced.
- **papers** contains relevant papers.
- **results** contains final figures and results.
- **src** organizes project analysis and preprocessing files into two main directories:
    - **gng:** Contains source files for the Go-NoGo task.
    - **rs:** Contains source files for the resting-state data, further subdivided into 20-second analyses and full recording analyses, and **notebooks,** which contains Jupyter notebooks for exploratory analysis.
- **work** contains Jupyter and Python files for exploratory analysis.
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
│   └── runs
├── docs
├── figures
├── papers
├── results
├── src
│   ├── gng
│   └── rs
│       ├── 20s
│       ├── full
│       │   ├── analysis
│       │   └── preprocessing
│       └── notebooks
└── work
```