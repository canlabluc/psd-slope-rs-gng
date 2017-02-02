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
- Code: BESA output, 11 means blink
- TriNo: Trigger Number. These are our port codes.

# Go-NoGo
This is the first set that we'll put through. Relevant trigger numbers and BESA codes are:
```
TrigNo
21 - Go prompt
22 - NoGo prompt
11 - Fixation
1  - Response (but BESA messed this up)

Code
11 - Blinks
```

