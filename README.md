# Directory and meta-data

This repository can be devided into four main parts: Experimental code (00), behavioural data (01), behavioural data analysis (02) and fMRI data analysis (03) for the "Mental Imagery of Needle-Stimulation Elicits Neural Activation Patterns with Acupuncture Point-Specificity" or short IMACU project.

## Directory:

### <b> 00_fMRI_experiment: </b> 
contains MATLAB code for running the experiment. Requirements: MATLAB 2018b or later, Psychtoolbox 3.0 or later
  - A_sound_cues.m: creates and saves the cues given to the acupuncturist during the fMRI experiment. Only needs to be run once
  - B_randomizer.m: balances cue shapes across participants, randomizes point labels across participants (thus used for blinding), randomizes colour assignment to points across runs. Needs to be run once per participant prior to the experiment. Creates and saves variables (save name: cue-assignment_subject-<sj_id>) used by all following scripts. These need to be loaded or kept in working memory prior to running any of the following.
  - C_point_test.m: promts the acupuncturist to stimulate the needles in the 3 points in order of how they were labeled for the respective participant (labels "1", "2", "3" where any of these might be the acupuncture point or control points, see randomizer) while informing the participant per visual text on the screen which point is being stimulated. After every point whas stimulates twice, the text instead promts the participant to answer out loud (experimenters recorded their answers) which point was stimulated. This was used to ensure that participants could discern the feeling per point.
  - D_test_run.m: a shorter version of E_experimental_script that did not require and MRI trigger and was run during the anatomical scan to get participants used to the task.
  - E_experimental_script.m: as described in the manuscript, saves log-files to be used during analysis

### <b> 01_behavioural_data: </b> 
contains log-files and MASS questionnaires as acquired before and after each session
  - A_Demographics: Identification, age (in years), weight (in kg), height (in cm) and sex (F = female) at time of acquisition
  - B_point_label_assignment: how the acupuncture point (Acu) and the control points (C1, C2) were labeled per participant
  - C_log_files: .tsv per participant, session and run, date-string marks time of acquisition. One row per trial with the following columns: TrialNr (condition, 1-3 for stimulation trials with number corresponding to participant-specific point labels, 4-6 for imagery trials),	onset (onset of the trial since first fMRI trigger in seconds),	preITI (inter-stimulus interval between stimulation/imagery and promt for response in seconds), Response (participants' answer of 1 to 4 per button press, 0 codes a missed response),	postITI (ITI until next trial in seconds)
  - D_MASS: Massachusetts Acupuncture Sensation Scale, modified to assess expectation of acupuncture sensations and expectation of sensations during imagery of acupuncture and sensations during acupuncture and imagery thereof. The MASS questionnaire has 12 items to rate on an anchored scale from 0 (none) to 10 (unbearable). The 13th is an open text field to fill in a non-listed sensation and its rating. The MASS questionnaires were administered in both sessions of measurement, before and after the experiment. Before the measurement, the question for participants was to rate their expectations about stimulation or imagery of acupuncture, for each point on the leg. Here, the questionnaires were preceded by the prompt "Please rate your expectations of today’s sensations during the acupuncture / imagery of acupuncture. Take a moment to evaluate what sensations you expect of todays session. If you expect something that is not listed, please write it on the blank marked “Other”."
After the measurement, they rated the sensations during stimulation or imagery of acupuncture, per each point. This was prompted by the statement "Please rate the sensations during today’s acupuncture stimulation / imagery of acupuncture stimulation. Take a moment to recall any sensations you experienced. If you felt something that is not listed, please write it on the blank marked “Other”."
On the first session both stimulation and imagery conditions were rated, and on the second session only imagery conditions were examined.

### <b> 02_behavioural_analysis: </b> 
contains python code to analyse behavioural data described above. Requirements: python3, numpy, pandas, matplotlib, scipy, os, statsmodels, jupyter notbook (recomended)
  - A_behav.ipynb/.py: Perfoms behavioural analyses (ANOVAs, etc) on data from log-files and MASS questionnaires
  - B_behav_figures.ipynb/.py: plots behavioural data as shown in mauscript

### <b> 03_fMRI_analysis: </b> 
all MATLAB code used to analyse the fMRI data. Can be devided into conversion (A), preprocessing (B), univariate analyses (C), multivariate analyses (D) and group inference (D). All require MATLAB 2023b or later
  - <b> A_conversion: </b> code used to convert dicom files to nifti files and sort them into BIDS (brain imaging data structure). Requirements: dicom2nifti toolbox ( https://de.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii )
    - A0_dicom2bids_batch.m: sets paths and requirements for conversion
    - A1_dicm2bids.m: executes conversion
  - <b> B_preprocessing: </b> code used to preprocess data. For univariate analyses, all preopocessing steps were run. For multivariate analyses only realignment and component correction in subject space. Requirements: SPM12 ( http://www.fil.ion.ucl.ac.uk/spm/ ), hMRI (https://www.cbs.mpg.de/departments/neurophysics/software/hmri-toolbox), for CompCorr: dPABI toolbox ( http://rfmri.org/dpabi )
    - B0_preprocessing_batch.m: batch administering the preprocessing steps per data set in specified order
    - B1_Realignment_all_runs.m: motion correction or realignment over all specified runs together. Saves realignment parameters as .txt to directory containing the nifits
    - B2_coregister_est.m: coregistration of anatomical scan onto functional scan, only changing meta-data (origin)
    - B3_segmentation.m: segmentation of anatomical image
    - B4_normalization.m: "warping" functional data based on anatomical image to MNI space
    - B5_smoothing.m: smoothing functional data with gaussian kernel of specified FWHM
    - B6_compcorr.m: component correction extracting a specified number of principal components from white matter and cerebrospinal fluid (masks in mni space or per participant in subject-space required). Saves principal components as .txt to directory containing the nifits
  - <b> C_univariate: </b>
