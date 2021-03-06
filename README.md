# Broken-Arm-Study

Processes imaging, accelerometry, and strength data from Broken Arm Study.

Key scripts:
- cast_proc.csh
  - Performs standard MRI preprocessing of T1, T2 and fMRI scans
  - Inputs:
    - DICOMS - in CNDA format
    - ${subject}_sessions.txt - name of sessions in CNDA download
    - ${subject}_sessions.txt - name of output sessions
    - instructions.txt - contains sequence parameters and processing flags
  - Outputs:
    - T1 - debiased, averaged across acquisitions, in native + atlas space
    - T2 - debiased, averaged across acquisitions, in native + atlas space
    - Freesurfer outputs - surfaces in native + atlas space, 128 + 32k
    - Preprocessed fMRI timecourses - debanded, mode 1000 normalized,
      slice-timing corrected, motion corrected, in atlas space, 4dfp format
    - RSFC Processed fMRI timecourses - bandbassed, scrubbed, nuisance
      regressed, in atlas space, 4dfp format
    - fMRI CIFTIs - Preprocessed (task) + RSFC processed (rest), 2D cortical
      surface + masked 3D subcortical timecourses
  - Must set download + output paths at top of script
  - Can run specific sections by uncommenting goto commands at top of script
- cast_analysis.m
  - Performs all secondary analyses
    - ANOVA-based event-related task analysis
    - Task-based ROI generation
    - ROI-ROI and ROI seed-based RSFC calculation
    - Longitudinal analysis of ROI-ROI RSFC measurements
    - Pulse detection + whole-brain pulse mapping
    - Accelerometry + strength testing
  - All inputs generated by cast_proc.csh

Dependencies (linked in src/)
  - 4dfp tools
  - FSL
  - Freesurfer
  - HCP workbench
  - Matlab
  - FieldTrip Matlab toolbox: ft_read_cifti_mod.m

How to run:
  - Edit instructions.txt to update sequence parameters + processing flags
  - Edit cast_proc.csh to update download and output paths
  - Run cast_proc.csh
  - Edit cast_analysis.m to update input and output paths
  - Run cast_analysis.m
