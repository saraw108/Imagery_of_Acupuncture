
%% decoding BATCH %%

% required toolboxes:
% hmri, spm

%#####################################################
%#################### INPUT ##########################
%#####################################################

clear runs


%SPM-path
SPM_path  = '.../toolboxes/spm12';

%data source directory
src_dir      = '.../IMACU/Data';

addpath(genpath('.../IMACU/C_statistical_analysis-main')); %%%%%%% your current working directory %%%%%%%%%%%%%% eg: 'F:\GForce\Dataanalysis\D_decoding_BIDS-main'
addpath(genpath('.../toolboxes/hMRI-toolbox'));
addpath('.../toolboxes/spm12');

%subject identifiers if all subjects are to be included
%%%% to do: SJs to sub
cd(src_dir)
clear SJs
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

display('Subjects found:')
SJs % analysis for these subjects
excludeSJ = [11]; % exclude those sjs:
ex_sub = 2;
excludeRuns = [];
expected_runs = 7;

jsons = 1;

%% selection of analysis steps (1-5) to be performed
analysis_switch = [ 1 2 ]; % 1 2 3 4 5
start_prefix='s8wr'; %eg. s8wra % s8wSEBsite1rb

zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files)
    for z = 1:size(zip_files, 1)
        gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
    end
end

%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7 %################################### skript can't (yet) process multiple sessions per subject #################################
    cd([src_dir filesep SJs{1}])
    sd = dir('ses*');
    sessNum = length(sd);
    for sess = 1:sessNum
        sessions(1, sess) = {sd(sess).name};
        for sb = 1:numel(SJs)
            cd([src_dir filesep SJs{sb} filesep sessions{sess} filesep 'func']);
            rd = dir('sub*.nii');
            for r = 1:length(rd)
                runs(sb, r, sess) = {rd(r).name};
            end
        end
    end
else
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep 'func']);
        rd = dir('sub*.nii');
        for r = 1:length(rd)
            runs(sb, r) = {rd(r).name};
        end
    end
end

%anatomy identifier
ana=['anat'];

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold.nii'])); % 55555555555555555555555555555555555555555555555
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%now we get the data from the json file
json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold.json']))); %extract all json files, althoguh they should have the same info
%because for some datasets (flicker and Ganzfeld) we have json files for
%each nift file and these are named differently, we have to check if the
%first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject
if isequal(size(json_files), [0, 1])
    json_files = (dir(fullfile(src_dir, '**', example_runs)));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order
slice_order = y';

%%% I suggest using the json values at least for TR since we know it is
%%% missing in the nifti headers for some datasets
n_slices = n_slices_json; % number of slices
TR=TR_json; % repetition time in sec.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1: extract onsets reeeealy fast if you havn't done that yet
% condition names
logDir='.../IMACU/Logs'; % in case you keep them somewhere differnt
log_pref = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  1st level glm for normalized data
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec

per_sess = 1;

include = [1]; % sessions to include

beta_dir = ['1st_level_BASIC_wr'];

condnames = {'StimAcu', 'StimC1', 'StimC2', 'ImagAcu', 'ImagC1', 'ImagC2', 'buttonPress'}; % {'ImagAcu', 'ImagC1', 'ImagC2', 'buttonPress'}; %
durations = 5; % epoch duration; for single events set 0
tr = TR;
fmri_t = 16; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0	= 1;
hpf      = 128; % High-pass filter cut-off; default 128 sec
% include multiple regressors (1=yes)
% if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st
hm=1;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=1;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)

mask = ''; % '' if implicit mask, otherwise specify

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 3:  contrasts
analysisfolder = [beta_dir '_hmcc']; % '_hm'

is_2nd = 0;

c_type = 't'; % F or t
cnames = {'StimAcu', 'StimC1', 'StimC2', 'ImagAcu', 'ImagC1', 'ImagC2'}; %{'ImagAcu', 'ImagC1', 'ImagC2'}; %

%very basic, adjust yourself, maybe model response regressor?
cvecs  = {[ 1  0  0  0  0  0  0 ], ...   % 1
    [ 0  1  0  0  0  0  0 ], ...   % 2
    [ 0  0  1  0  0  0  0 ], ...   % 3
    [ 0  0  0  1  0  0  0 ], ...   % 4
    [ 0  0  0  0  1  0  0 ], ...   % 5
    [ 0  0  0  0  0  1  0 ]};      % 6

del=1; % Delete existing contrasts (1=yes)
n_hm = 6;
n_cc = 5;

covar = 0;

%%
currPrefix = start_prefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        %% extract onsets
        case 1

            for s = 1:numel(SJs)

                if ismember(s, excludeSJ)
                    continue;
                else
                    cd(logDir)

                    for ses = 1:sessNum
                        oruns=dir([log_pref 'LogFile_subject_' num2str(s) '_sess_' num2str(ses) '_*.tsv']); %do they have names??

                        if ~isempty(runs)

                            for r=1:size(runs, 2)

                                this_log = tdfread([logDir filesep oruns(r).name]);

                                % cond 1: Stimulation on point Acu
                                onsets{s,ses,r,1} = this_log.onset(this_log.TrialNr == 1) + 1; %%% change to struct with onsets.irgendwas{
                                % cond 2: stim C1
                                onsets{s,ses,r,2} = this_log.onset(this_log.TrialNr == 2) + 1;
                                % cond 3: stim C2
                                onsets{s,ses,r,3} = this_log.onset(this_log.TrialNr == 3) + 1;
                                % cond 4: Imag Acu
                                onsets{s,ses,r,4} = this_log.onset(this_log.TrialNr == 4) + 1;
                                % cond 5: Imag C1
                                onsets{s,ses,r,5} = this_log.onset(this_log.TrialNr == 5) + 1;
                                % cond 6: Imag C2
                                onsets{s,ses,r,6} = this_log.onset(this_log.TrialNr == 6) + 1;

                                % regr 7: Stim resp
                                onsets{s,ses,r,7} = this_log.onset + 1 + this_log.preITI;
                                % regr 8: Imag resp


                                % cond 1: Stimulation on point Acu
                                responses(s,ses,r,1) = 1; % mean(this_log.Response(this_log.TrialNr == 1)); %%% change to struct with onsets.irgendwas{
                                % cond 2: stim C1
                                responses(s,ses,r,2) = 1; % mean(this_log.Response(this_log.TrialNr == 2));
                                % cond 3: stim C2
                                responses(s,ses,r,3) = 1; % mean(this_log.Response(this_log.TrialNr == 3));
                                % cond 4: Imag Acu
                                responses(s,ses,r,4) = 1; % mean(this_log.Response(this_log.TrialNr == 4));
                                % cond 5: Imag C1
                                responses(s,ses,r,5) = 1; % mean(this_log.Response(this_log.TrialNr == 5));
                                % cond 6: Imag C2
                                responses(s,ses,r,6) = 1; % mean(this_log.Response(this_log.TrialNr == 6));

                                % regr 7: Stim resp
                                responses(s,ses,r,7) = 1;
                                % regr 8: Imag  resp

                            end
                        end %end if
                    end
                end
            end

            %% GLM
        case 2

            if sessNum == 0
                % folder that will contain the created job.mat file and SPM file

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        display(['Step 2, 1st level glm: ' SJs{sj} ])
                        subj_dir = fullfile(src_dir, SJs{sj});
                        these_onsets = reshape(onsets(sj,1,:,:), size(runs,2), length(condnames));

                        D2_1st_level_glm(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, these_onsets, durations, hm, cc, mask);
                    end
                end
            elseif per_sess == 1

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        for ses = include

                            display(['Step 2, 1st level glm: ' SJs{sj} ])
                            ses_dir = fullfile([src_dir filesep SJs{sj}], ['ses-' num2str(ses)]);
                            these_onsets = reshape(onsets(sj,ses,:,8-length(condnames):7), size(runs,2), length(condnames));
                            these_responses = reshape(responses(sj,ses,:,8-length(condnames):7), size(runs,2), length(condnames));
                            if tr~= 1
                                tr = 1;
                            end
                            D2_1st_level_glm(SJs, sj, ses_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, reshape(runs(:,:,ses), numel(SJs), size(runs,2)), condnames, these_onsets, durations, hm, cc, mask, these_responses);

                        end
                    end
                end
            else

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        display(['Step 2, 1st level glm: ' SJs{sj} ])
                        subj_dir = fullfile(src_dir, SJs{sj});
                        these_onsets = reshape(onsets(sj,:,:,:), size(runs,2)*2, length(condnames));
                        these_responses = reshape(responses(sj,ses,:,:), size(runs,2), length(condnames));
                        covar = these_responses./sum(these_responses,1)

                        ses_dir = spm_select('List', subj_dir, 'dir', ['ses-[' erase(num2str(include), ' ') ']']);
                        D2_1st_level_glm(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, durations, ses_dir, hm, cc, mask, covar);
                    end
                end
            end
    		%% contrasts
        case 3
            if sessNum == 0
                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        display(['Step 3, 1st level contrasts: ' SJs{sj} ])
                        subj_dir = fullfile(src_dir, SJs{sj});
                        C3_contrast_1stLevel(subj_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs, is_2nd, c_type);
                    end
                end
            elseif per_sess == 1

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        for ses = include

                            display(['Step 3, 1st level contrasts: ' SJs{sj} ])
                            ses_dir = fullfile([src_dir filesep SJs{sj}], ['ses-' num2str(ses)]);

                            if covar == 0
                                C3_contrast_1stLevel(ses_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs(:,:,ses), is_2nd, c_type);
                            else
                                these_responses = reshape(responses(sj,ses,:,4:7), size(runs,2), length(condnames));
                                covar = these_responses./sum(these_responses,1);

                                C3_contrast_1stLevel_covar(ses_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs(:,:,ses), is_2nd, c_type, covar);
                            end
                        end
                    end
                end
            else

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        display(['Step 3, 1st level glm: ' SJs{sj} ])
                        subj_dir = fullfile(src_dir, SJs{sj});
                        these_onsets = reshape(onsets(sj,:,:,:), size(runs,2)*2, length(condnames));

                        ses_dir = spm_select('List', subj_dir, 'dir', ['ses-[' erase(num2str(include), ' ') ']']);
                        C3_contrast_1stLevel(subj_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs, is_2nd, c_type);
                    end
                end
            end
    end
end
disp('Done.')

