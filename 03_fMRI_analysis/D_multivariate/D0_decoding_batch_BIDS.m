%% decoding BATCH %%

% required toolboxes:
% hmri, spm
% the decoding toolbox TDT
% https://doi.org/10.3389/fninf.2014.00088

%#####################################################
%#################### INPUT ##########################
%#####################################################


%SPM-path
SPM_path  = '.../toolboxes/spm12';

%data source directory
src_dir      = '.../IMACU/Data';

addpath(genpath('.../IMACU/decoding-main'));
addpath(genpath('.../toolboxes/hMRI-toolbox'));
addpath('.../toolboxes/decoding_toolbox');
addpath('../toolboxes/spm12');

%subject identifiers if all subjects are to be included
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

display('Subjects found:')
SJs % analysis for these subjects
display('Subjects to keep:')
excludeSJ = [11]; % delete those sjs

zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files)
    for z = 1:size(zip_files, 1)
        gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
    end
end

%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7
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

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%now we get the data from the json file
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info
if isequal(size(json_files), [0, 1])
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold.json'])));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order
slice_order = y';

%now get the same info from nifti header
nifti_file_metadata = [nifti_files(1).folder, filesep, nifti_files(1).name];
info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4);
n_slices_nifti = info.ImageSize(3);
vox_size=repmat(info.PixelDimensions(1),1,3);

%compare json and nifti header
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti")
end
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti")
end

n_slices = n_slices_json; % number of slices
TR=TR_json; % repetition time in sec.

%% selection of analysis steps (1-5) to be performed
analysis_switch = [ 1 2 3 4 5 6 7 ]; % 1 2 3 4 5
start_prefix='r';
%'' when bids and unpreprocessed;
%'r' when realigned
%'' if doing normilazation of accuracy maps
% if aready normalized
% watch order of prefixes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1: extract onsets reeeealy fast if you havn't done that yet
logDir='.../IMACU/Logs'; % in case you keep them somewhere differnt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2: 1st level glm for NOT-normalized data to decode on!
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec

per_sess = 1; % for each session separately? 1=yes
include = [1]; % sessions to include

beta_dir = ['1st_level_DEC_r']; % I recommend adding the prepro prefix to avoid confusion with univariate analyses :)

condnames = {'StimAcu', 'StimC1', 'StimC2'; 'ImagAcu', 'ImagC1', 'ImagC2'}; % {'ImagAcu', 'ImagC1', 'ImagC2'}; %
durations = 5; % stim duration; for single events set 0
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
%# step 3:  decode (regular)

over_sess = 0;

map_type_4 = {'accuracy_minus_chance', 'confusion_matrix'}; % results to write, options:
% 'accuracy', 'accuracy_minus_chance', 'sensitivity_minus_chance', 'specificity_minus_chance', 'balanced_accuracy_minus_chance', 'confusion_matrix', 'AUC_minus_chance'

from_to_dec = [1 2]; % only used if over_sess = 1, from which session to which? e.g. [1 2]
include_dec = [1 2];

betas = ['1st_level_DEC_r_hmcc'];
grid_C = 0; % parameter optimization? 0 = no

dec_type = "searchlight"; % searchlight or roi
dec_task = "classification"; % regression (SVR) and classification (SVM) are implemented so far

labelnames 	= 	{'StimAcu', 'StimC1', 'StimC2'; 'ImagAcu', 'ImagC1', 'ImagC2'};  % 'ImagAcu', 'ImagC1', 'ImagC2'
labels = [1 2 3];

Nbins = 1; % default 1
con_pairs = [1 2];

sbj_level_folder = 'DEC';
outputDir4 = 'D'; % folder prefix

rad = 3; % radius of searchlight

mask_ident = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 4: xclass

over_ses = 0;

map_type_X = {'accuracy_minus_chance', 'confusion_matrix'}; %

to_include_x = [1];

from_to = [2 1]; % only used if over sess = 1: e.g. decode from session 1 to sess 2 [1 2] or other way around [2 1]

outputDirX = 'XClass';
% Set the label name pairess to the regressor names which you want to use for decoding, e.g. 'higher' and 'lower'
% if you don't remember the names --> run display_regressor_names(betaDir)
% labels of training and test data must match for proper cross classification!

labelnames_train = {'StimAcu', 'StimC1', 'StimC2'; 'ImagAcu', 'ImagC1', 'ImagC2'};
labelnames_test  = {'ImagAcu', 'ImagC1', 'ImagC2'; 'StimAcu', 'StimC1', 'StimC2'};

radX = 3; % radius of searchlight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 5: rewrite confusion to 3D

cell_titles = {'TrueAcu',	'AcuToC1',	'AcuToC2',	'C1ToAcu',	'TrueCp1',	'Cp1ToC2',	'C2ToAcu',	'Cp2ToC1',	'TrueCp2'};
decoding_folders = {'ses-1/DEC/D_StimAcu_StimC1_StimC2', 'ses-1/DEC/D_ImagAcu_ImagC1_ImagC2', 'ses-2/DEC/D_ImagAcu_ImagC1_ImagC2', 'ses-1/DEC/XClass_StimAcu_to_ImagAcu', 'ses-1/DEC/XClass_ImagAcu_to_StimAcu'}; % or 'DEC/D_ImagAcu_ImagC1_ImagC2_ses_1_to_2' or whatever

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 6:  normalize

map_type = 'percent_.*\'; % 'percent_.*\' in case of confusion matrix results % 'res_accuracy_minus_chance' otherwise

include_norm_smooth = [1]; % sessions to perform this on

over_sess_norm = 1; % was realignment over sessions?
over_ses_dec = 0; % was decoding over sessions?

folder_inst1 = 'DEC';
folder_inst2 = {'D_ImagAcu_ImagC1_ImagC2'}; % or ... 'D_StimAcu_StimC1_StimC2',  'D_ImagAcu_ImagC1_ImagC2_ses_1_to_2', 'D_ImagAcu_ImagC1_ImagC2_ses_2_to_1'
% XClass_ImagAcu_vs_ImagC1_to_StimAcu_vs_StimC1, 'XClass_StimAcu_to_ImagAcu', 'XClass_ImagAcu_to_StimAcu'

vox_size = [2 2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 7:  smooth
s_kernel = [5 5 5];
% takes input from above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currPrefix=start_prefix;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        case 1

            for s = 1:numel(SJs)

                if ismember(s, excludeSJ)
                    continue;
                else
                    % subj_dir = fullfile(src_dir, SJs{s});
                    cd(logDir)

                    for ses = 1:sessNum
                        oruns=dir(['LogFile_subject_' num2str(s) '_sess_' num2str(ses) '_*.tsv']); %do they have names??

                        if ~isempty(runs)

                            for r=1:size(runs, 2)

                                this_log = tdfread([logDir filesep oruns(r).name]);

                                % cond 1: Stimulation on point Acu
                                onsets{s,ses,r,1} = this_log.onset(this_log.TrialNr == 1) + 1;
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

                            end
                        end
                    end
                end
            end

        	%% glm
        case 2

            if sessNum == 0

                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        display(['Step 3, 1st level glm: ' SJs{sj} ])
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

                            display(['Step 3, 1st level glm: ' SJs{sj} ])
                            ses_dir = fullfile([src_dir filesep SJs{sj}], ['ses-' num2str(ses)]);
                            these_onsets = reshape(onsets(sj,ses,:,(7-length(condnames)):6), size(runs,2), length(condnames));

                            D2_1st_level_glm(SJs, sj, ses_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, reshape(runs(:,:,ses), numel(SJs), size(runs,2)), condnames, these_onsets, durations, hm, cc, mask);

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
                        D2b_1st_level_glm(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, durations, ses_dir, hm, cc, mask);

                    end
                end
            end


        	%% decoding
        case 3

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    if over_sess == 1

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            masks = 1;
                        end


                        SJ_dir = [src_dir filesep SJs{sj}];
                        dec_folder = [SJ_dir filesep sbj_level_folder];
                        ses_dir = spm_select('List', SJ_dir, 'dir', ['ses-[' erase(num2str(include_dec), ' ') ']']);

                        if ~exist(dec_folder, 'dir')
                            mkdir(dec_folder);
                        end

                        if dec_task == "classification"

                            for cp = 1:size(labelnames,1)

                                label1 = labelnames(cp,1);
                                label2 = labelnames(cp,2);
                    			label3 = labelnames(cp,3);
                                these_labelnames = [label1, label2, label3];

                                display(['Step 3, Decoding: ' SJs{sj} ', Conditions: ' cell2mat(label1) ', ' cell2mat(label2) ', ' cell2mat(label3)]) %
                                beta_path    = strcat(ses_dir, [filesep betas]);
                                beta_paths    = strcat(SJ_dir, filesep, beta_path)
                                output_path  = fullfile(dec_folder, [outputDir4 '_' cell2mat(label1) '_' cell2mat(label2) '_' cell2mat(label3) '_ses_' num2str(from_to_dec(1)) '_to_' num2str(from_to_dec(2))]); %

                                D3b_Decoding_batch_xSess(beta_paths, output_path, these_labelnames, these_labelnames, from_to_dec, dec_type, rad, masks, map_type_4, grid_C);

                                clear label1 label2 label3 these_labelnames

                            end
                        end

                    elseif sessNum > 0 % but per session a thing

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            masks = 1;
                        end

                        for r = 1:length(masks)

                            SJ_dir = [src_dir filesep SJs{sj}];
                            for ses = include_dec

                                ses_path = [SJ_dir filesep 'ses-' num2str(ses)]; %%%%%%
                                dec_folder = [ses_path filesep sbj_level_folder];
                                if ~exist(dec_folder, 'dir')
                                    mkdir(dec_folder);
                                end

                                for b = 1:Nbins

                                    if dec_task == "classification"

                                        for cp = 1:size(labelnames,1) % maybe adjust to allow restriction of used betas beforehand

                                            label1 = labelnames(cp,1);
                                            label2 = labelnames(cp,2);
                            			    label3 = labelnames(cp,3);
                                            these_labelnames = [label1, label2, label3];

                                            display(['Step 3, Decoding: ' SJs{sj} ', Conditions: ' cell2mat(label1) ', ' cell2mat(label2) ', ' cell2mat(label3)]) %
                                            beta_path    = fullfile(ses_path, betas);
                                            output_path  = fullfile(dec_folder, [outputDir4 '_' cell2mat(label1) '_' cell2mat(label2) '_' cell2mat(label3)]); %
                                            D3class_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks, map_type_4, labels, grid_C);

                                            clear label1 label2 these_labelnames

                                        end

                                    end

                                end
                            end
                        end

                    else % usual

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            cd(src_dir);
                            masks = 1;
                        end

                        for r = 1:length(masks)

                            SJ_dir = [src_dir filesep SJs{sj}];
                            dec_folder = [SJ_dir filesep sbj_level_folder];
                            if ~exist(dec_folder, 'dir')
                                mkdir(dec_folder);
                            end


                            for b = 1:Nbins

                                if dec_task == "classification"

                                    for cp = 1:size(labelnames,1) % maybe adjust to allow restriction of used betas beforehand

                                        label1 = labelnames(cp,1);
                                        label2 = labelnames(cp,2);
                            			label3 = labelnames(cp,3);
                                        these_labelnames = [label1, label2, label3];

                                        display(['Step 3, Decoding: ' SJs{sj} ', Conditions: ' cell2mat(label1) ', ' cell2mat(label2) ', ' cell2mat(label3)]) %
                                        beta_path    = fullfile(betas);
                                        output_path  = fullfile(dec_folder, [outputDir4 '_' cell2mat(label1) '_' cell2mat(label2)]);
                                        D3class_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks, map_type_4, grid_C);
                                        clear label1 label2 label3 these_labelnames

                                    end

                                end

                            end
                        end

                    end

                end

            end


        %% decoding XClass
        case 4

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    if over_ses == 1

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            cd(src_dir);
                            masks = 1;
                        end

                        SJ_dir = [src_dir filesep SJs{sj}];

                        % ses_path = [SJ_dir filesep 'ses-' num2str(ses)]; %%%%%%
                        ses_dir = spm_select('List', SJ_dir, 'dir', ['ses-[' erase(num2str(to_include_x), ' ') ']']);
                        dec_folder = [SJ_dir filesep sbj_level_folder];
                        if ~exist(dec_folder, 'dir')
                            mkdir(dec_folder);
                        end

                        for b = 1:Nbins

                		    for cp = 1:size(labelnames_train,1) % maybe adjust to allow restriction of used betas beforehand

                                these_train_labels = labelnames_train(cp, :);
                                these_test_labels = labelnames_test(cp, :);

                                display(['Step 4, xDecoding: ' SJs{sj} ', Conditions: ' outputDirX 'from ses ' num2str(from_to(1)) ' ' cell2mat(these_train_labels(1)) ' to ses ' num2str(from_to(2)) ' ' cell2mat(these_test_labels(1))])
                                beta_path    = strcat(ses_dir, [filesep betas]);
                                beta_paths    = strcat(SJ_dir, filesep, beta_path);
                                output_path  = fullfile(dec_folder, [outputDirX '_' cell2mat(these_train_labels(1)) '_to_' cell2mat(these_test_labels(1)) '_ses_' num2str(from_to(1)) '_to_' num2str(from_to(2))]);

                                D4b_Decoding_batch_xClass(beta_paths, output_path, these_train_labels, these_test_labels, from_to, dec_type, radX, masks, map_type_X, grid_C)% hier is das andere wichtige

                                clear label1 label2 these_labelnames
                		    end

                        end

                    elseif sessNum > 0 % but per session a thing

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            cd(src_dir);
                            masks = 1;
                        end

                        SJ_dir = [src_dir filesep SJs{sj}];
                        for ses = to_include_x

                            ses_path = [SJ_dir filesep 'ses-' num2str(ses)]; %%%%%%
                            dec_folder = [ses_path filesep sbj_level_folder];
                            if ~exist(dec_folder, 'dir')
                                mkdir(dec_folder);
                            end

                            for cp = 1:size(labelnames_train,1) % maybe adjust to allow restriction of used betas beforehand

                                these_train_labels = labelnames_train(cp, :);
                                these_test_labels = labelnames_test(cp, :);

                                display(['Step 4, xDecoding: ' SJs{sj} ', Conditions: ' cell2mat(these_train_labels(1)) ' to ' cell2mat(these_test_labels(1))])
                                beta_path    = fullfile(ses_path, betas);
                                output_path  = fullfile(dec_folder, [outputDirX '_' cell2mat(these_train_labels(1)) '_to_' cell2mat(these_test_labels(1))]);

                                D4_Decoding_batch_xClass(beta_path, output_path, these_train_labels, these_test_labels, dec_type, radX, masks, map_type_X, grid_C)

                            end
                        end

                    else % usual

                        if dec_type == "roi"
                            masks = dir(['^' mask_ident '*.nii']); %or whatever you use as an identifier
                        else
                            cd(src_dir);
                            masks = 1;
                        end

                        for r = 1:length(masks)

                            SJ_dir = [src_dir filesep SJs{sj}];
                            for ses = 1:sessNum

                                sub_path = SJ_dir; %%%%%%
                                dec_folder = [ses_path filesep sbj_level_folder num2str(r)];
                                if ~exist(dec_folder, 'dir')
                                    mkdir(dec_folder);
                                end

                                for cp = 1:size(labelnames_train,1) % maybe adjust to allow restriction of used betas beforehand

                                    these_train_labels = labelnames_train(cp, :);
                                    these_test_labels = labelnames_test(cp, :);

                                    display(['Step 4, xDecoding: ' SJs{sj} ', Conditions: ' cell2mat(these_train_labels{1}) ', ' cell2mat(these_test_labels{1})])
                                    beta_path    = fullfile(sub_path, betaDir);
                                    output_path  = fullfile(dec_folder, [outputDirX '_' cell2mat(these_train_labels(1)) '_' labelnames_test{2}]);

                				    D4_Decoding_batch_xClass(beta_path, output_path, these_train_labels, these_test_labels, dec_type, radX, masks, map_type_X)

                                end
                            end
                        end
                    end
                end

            end

    	case 5 %% calculate 3D confusion

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    sj_dir = [src_dir filesep SJs{sj}];

                    for fold = decoding_folders

                        working_dir = [sj_dir filesep fold{:}];

                        D5_confusion_to3D(working_dir, 'res_confusion_matrix.mat', 'res_accuracy_minus_chance.nii', cell_titles);

                    end

                end

            end

        case 6 %% norm

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    if sessNum > 0

                        if over_sess_norm == 1 % if realignment was done over sessions

                            sj_dir = [src_dir filesep SJs{sj}];
                            ses_dir = spm_select('List', sj_dir, 'dir', ['ses-[' erase(num2str(1:sessNum), ' ') ']']);
                            struct_dir = [sj_dir filesep ses_dir(1,:) filesep 'anat'];
                            func_dir = [sj_dir filesep ses_dir(1,:) filesep 'func'];

                            Images = [];


                            if over_ses_dec == 1 % if decoding results are on sj instance

                                for fold2 = folder_inst2

                                    data_dir = [sj_dir filesep folder_inst1 filesep fold2{:}];

                					f3 = spm_select('List', data_dir, ['^' map_type '.nii']);
                                    numVols = size(f3,1);
                                    Images{1} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);

                                    D6c_normalization(Images, struct_dir, vox_size);
                                end

                            else % if decoding results are on session instance

                                for ses = include_norm_smooth

                                    sj_dir = [src_dir filesep SJs{sj}];
                                    this_dir = [sj_dir filesep ses_dir(ses,:)];

                                    for fold2 = folder_inst2
                                        data_dir = [this_dir filesep folder_inst1 filesep fold2{:}];
                                        f3 = spm_select('List', data_dir, ['^' map_type '.nii']);
                                        if isempty(f3)
                                            f3 = spm_select('List', data_dir, [map_type '.nii']);
                                        end
                                        numVols = size(f3,1);
                                        Images{1} = cellstr([repmat([data_dir filesep], numVols, 1) char(f3) repmat(',1', numVols, 1)]);

                                        D6c_normalization(Images, struct_dir, vox_size);
                                    end

                                end

                            end

                        else % if realignment was NOT done over sessions but separatly

                            for ses = 1:sessNum

                                this_dir = [sj_dir filesep ses_dir(ses,:)];

                                func_dir = [this_dir filesep 'func'];
                                struct_dir = [this_dir filesep 'anat'];

                                D6a_coregister_est(func_dir, struct_dir, Images);
                                D6b_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');

                                for fold2 = folder_inst2
                                    data_dir = [this_dir filesep folder_inst1 filesep fold2{:}]; 

                                    f3 = spm_select('List', data_dir, ['^' map_type '.nii']);
                                    if isempty(f3)
                                        data_dir = [sj_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                                        f3 = spm_select('List', data_dir, '^res_zcorr.nii');
                                    end
                                    numVols = size(f3,1);
                                    Images{1} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);

                                    D6c_normalization(Images, struct_dir, vox_size);
                                end

                            end

                        end

                    else % no sessions

                        sj_dir = [src_dir filesep SJs{sj}];
                        func_dir = [sj_dir filesep 'func'];
                        struct_dir = [sj_dir filesep 'anat'];

                        for b = 1:Nbins

                            data_dir = [sj_dir filesep folder_inst1 filesep folder_inst2]; % C:\Users\saraw\Desktop\evenmoremodules\theSeminar\bicep\Data\sub-001\DEC\D4_ALT_SIM
                            f3 = spm_select('List', data_dir, ['^' map_type '.nii']);
                            if isempty(f3)
                                data_dir = [sj_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                                f3 = spm_select('List', data_dir, '^res_zcorr.nii');
                            end
                            numVols = size(f3,1);
                            Images{b} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);

                        end

                        D6a_coregister_est(func_dir, struct_dir, Images);
                        D6b_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                        D6c_normalization(Images, struct_dir, vox_size);

                    end

                end

            end

            currPrefix = 'w';

        case 7 %% smooth

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    sj_dir = [src_dir filesep SJs{sj}];

                    if sessNum > 0 && over_ses_dec == 0

                        for ses = 1:include_norm_smooth

                            this_dir = [sj_dir filesep ses_dir(ses,:)];

                            for fold2 = folder_inst2

                                data_dir = [this_dir filesep folder_inst1 filesep fold2{:}];
                                D7_smoothing(data_dir, SJs{sj}, ['^' currPrefix map_type '.nii'], s_kernel);

                            end

                        end

                    else

                        for b = 1:Nbins

                            for fold2 = folder_inst2

                                data_dir = [sj_dir filesep folder_inst1 filesep fold2{:}];
                                D7_smoothing(data_dir, SJs{sj}, ['^' currPrefix map_type '.nii'], s_kernel);

                            end

                        end
                    end

                end

            end

    end
end
disp('Done.')

