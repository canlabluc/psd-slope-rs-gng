# PSD Slopes: Resting state, Go-NoGo data
This project organizes efforts to compute the differences in power spectral density (PSD) slope between aMCIs, SAs, and their controls, for both task-related (Go-NoGo, referred to as `gng`) and resting-state data (referred to as `rs`). PSD slope has previously been identified as a proxy for the measurement of neural noise by [Voytek and colleagues](http://www.jneurosci.org/content/35/38/13257). The slope of  PSDs are computed at both the sensor and source level (source modeling is done through BESA).

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
- **results** contains final figures and results.
- **src** organizes project analysis and preprocessing files into two main directories:
    - **gng:** Contains source files for the Go-NoGo task.
    - **rs:** Contains source files for the resting-state data, further subdivided into 20-second analyses and full recording analyses.

## running an analysis
Instructions on estimating neural noise on a set during resting state are located in the docs:
- `docs/Resting state pipelind.md`: For estimating neural noise on resting state data.
- `docs/Go-NoGo pipeline.md`: For estimating neural noise on Go-NoGo data.

