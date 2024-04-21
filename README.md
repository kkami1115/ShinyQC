# ShinyQC
GUI that automatically pairs fastq files in a specified directory and executes Rfastp with paired-end.

# Usage.
`fastq_dir` ... The directory where fastqs exists.

`result_dir` ... The directory where results (html,json) and processed fastqs are output.

`file_pattern` ... Read1 and Read2 patterns in regular expression. Example: `_R[12]`

When all files have been processed, the contents of the Overview and Visualize per sample tabs can be displayed.


# ToDo
- Match the visualize per sample with the actual html (the specification is slightly different now)
- Make duplicate rate visible in Overview
- Add settings for arguments passed to `Rfastp` in the Settings tab
- add a stop button during processing
