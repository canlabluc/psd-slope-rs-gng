# BESA - Task-related Pipeline

This markdown document covers the pipeline for putting the task-related data through the spectral slope analysis. There are two sets of task-related data:
- Go-NoGo task
- Visual categorization task

# `.evt` file format
BESA outputs event files as `.evt`s. An example `.evt` file:
```
Tmu             Code    TriNo   Comnt
2633203         11      0       Pattern 1
5355859         11      0       Pattern 1
12109766        11      0       Pattern 1
20594141        11      0       Pattern 1
29332422        11      0       Pattern 1
33773828        11      0       Pattern 1
39998438        11      0       Pattern 1
...             ...     ...     ...
```

Where:
- Tmu: Event time in microseconds
- Code: BESA output, from off-line processing.
- TriNo: Trigger Number. These are our port codes, added to the recording during the task.
- Comnt: Relevant event comment. Not really necessary for us to consider.

# Go-NoGo
Relevant trigger numbers and BESA codes are:
```
TrigNo
21 - Go prompt
22 - NoGo prompt
11 - Fixation
1  - Response (but BESA messed this up)

Code (wiki.besa.de/index.php?title=Event_File_Format)
11 - Blinks
21 - Artifact start
22 - Artifact stop
```

