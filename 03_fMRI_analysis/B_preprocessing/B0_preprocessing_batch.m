%### step A, structure data before running the preprocessing
%### --> convert DICOM images to 4D NIFTI image
%### --> subject folders including functional (e.g. 'func') and anatomy (e.g. 'anat') folders
%### --> optional instance for sessions

% required toolboxes:
% SPM12 ( http://www.fil.ion.ucl.ac.uk/spm/ )
% hMRI (https://www.cbs.mpg.de/departments/neurophysics/software/hmri-toolbox)
% CompCorr only: dPABI toolbox ( http://rfmri.org/dpabi )

%#####################################################
%#################### INPUT ##########################
%#####################################################


addpath('../IMACU/task-based_fMRI_preprocessing-main')
addpath(genpath('.../toolboxes/hMRI-toolbox'))
addpath(genpath('.../toolboxes/DPABI-master'))
addpath('.../toolboxes/spm12')

%SPM-path
SPM_path  = '/home/saraw80/toolboxes/spm12';

%data source directory
src_dir      = ['.../IMACU/Data'];

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

excludeSJ = [11];

%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7
    cd([src_dir filesep SJs{1}])
    sd = dir('ses*');
    sessNum = length(sd);
    for sess = 1:sessNum
        sessions(1, sess) = {sd(sess).name};
    end
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep sessions{1} filesep 'func']);
        rd = dir('sub*.nii');
        for r = 1:length(rd)
            runs(sb, r) = {rd(r).name};
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
%unzip
zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files)
    for z = 1:size(zip_files, 1)
        gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
    end
end

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%anatomical masks
mni_wm_mask=['.../toolboxes/wm_mask_eroded.nii']; %white matter mask file
mni_csf_mask=['.../toolboxes/csf_mask_eroded.nii']; %csf mask file
mni_full_brain_mask=['.../toolboxes/full_brain_mask.nii']; %full brain mask file

% selection of analysis steps (1-6) to be performed
analysis_switch = [1 2 3 4 5 6]; % in case of decoding afterwards: just [1 6] with comp-corr in sj-space
start_prefix=''; 

%now we get the data from the json file
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files
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

%compare json and nifti header
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti")
end
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti")
end

TR = TR_json;
n_slices = n_slices_json;

%% descriptions and settings of steps
%# step 1: Realignment                                --> prefix: r
rb = 1;                                     %%% if 1, realign over all runs
over_sess = 0;

%# step 2:  Coregister (estimate) mean-epi 2 anatomy
corrPrefix = '';

%# step 3: segmentation

%# step 4:  Normalization                              --> prefix: w
vox_size=[2 2 2]; % voxel size in mm

%# step 10 Smoothing                                   --> prefix: s
kernel_size=[8 8 8]; %FWHM kernel size

%# step 9 Calculate WM and CSF Nuisance Signal
numComp = 5; % number of principle components
sj_space = 0; % mni space (0) or sj-space (1)? latter requires computation of masks per sj (see below)


%#####################################################
%#################### INPUT end ######################
%#####################################################

%%
currPrefix=start_prefix;

for n=analysis_switch

    switch n

        case 1 %% Realignment
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    sj_dir = [src_dir filesep SJs{sj}];
                    if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                        display(['Step 1, realignment: ' SJs{sj} ', ' runs{sj, r}])
                        funcPath = [src_dir filesep SJs{sj} filesep 'func'];
                        run_dir = fullfile(funcPath);
                        for r = 1:size(runs, 2)
                            run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                        end

                        B1_Realignment_all_runs(sj_dir, run_files, over_sess);

                    elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}]) && over_sess == 0
                        for ses = 1:sessNum
                            display(['Step 1, realignment: ' SJs{sj} ', ' runs{sj, r}])
                            sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                            run_dir = fullfile(sesPath);
                            for r = 1:size(runs, 2)
                                run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                            end

                            B1_Realignment_all_runs(sesPath, run_files, over_sess);

                        end
                    elseif over_sess == 1 
                        cd([src_dir filesep SJs{sj}]);
                        runsies = dir(['**' filesep currPrefix SJs{sj} '*bold.nii']);
                        B1_Realignment_all_runs([src_dir filesep SJs{sj}], runsies, over_sess);

                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['r' currPrefix]; 
            

        case 2 %% Coregister (estimate)
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    if exist([src_dir filesep SJs{sj} filesep 'func'])
                        display(['Step 2, coregistration: ' SJs{sj}])
                        funcPath = [src_dir filesep SJs{sj}];
                        func_dir        = fullfile(funcPath, 'func');
                        struct_dir      = fullfile(funcPath, ana);
                        B2_coregister_est(currPrefix, func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
                    elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                        if over_sess == 0
                            for ses = 1:sessNum
                                display(['Step 2, coregistration: ' SJs{sj}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                func_dir        = fullfile(sesPath, 'func');
                                struct_dir      = fullfile(sesPath, ana);
                                B2_coregister_est(currPrefix, func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
                            end
                        else
                            display(['Step 2, coregistration: ' SJs{sj}])
                            sesPath = [src_dir filesep SJs{sj} filesep 'ses-1'];
                            func_dir        = fullfile(sesPath, 'func');
                            struct_dir      = fullfile(sesPath, ana);
                            B2_coregister_est(currPrefix, func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
                        end

                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} 's functional data do not exsist ###########'])
                    end
                end
            end

        case 3 %% Segmentation
            warning off
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                elseif exist([src_dir filesep SJs{sj} filesep ana])==7
                    display(['Step 3, segmentation: ' SJs{sj}])
                    struct_dir = fullfile(src_dir, SJs{sj}, ana);

                    B3_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');

                elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep ana])==7
                    for ses = 1:sessNum
                        sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep ana];
                        display(['Step 3, segmentation: ' SJs{sj}])
                        struct_dir = fullfile(sesPath);

                        B3_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');

                    end
                else
                    display('###########################################################')
                    display(['############### ' SJs{sj} ', ' ana ' does not exsist ###########'])
                end
            end
        
        case 4 %% Normalization
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else

                    if exist([src_dir filesep SJs{sj} filesep 'func'])
                        display(['Step 4, normalization: ' SJs{sj}])
                        funcPath = [src_dir filesep SJs{sj}];
                        struct_dir = fullfile(funcPath, ana);
                        data_dir = fullfile(funcPath, 'func');
                        B4_normalization(data_dir, struct_dir, sj, runs, vox_size, currPrefix);
                    elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                        for ses = 1:sessNum
                            display(['Step 4, normalization: ' SJs{sj}])
                            sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                            func_dir        = fullfile(sesPath, 'func');
                            struct_dir      = fullfile(sesPath, ana);
                            B4_normalization(data_dir, struct_dir, sj, runs, vox_size, currPrefix);
                        end
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                    end
                end
            end
            currPrefix=['w' currPrefix];

        case 5 %% Smoothing
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 5, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            run_dir = fullfile(funcPath, 'func');
                            B5_smoothing(run_dir, SJs{sj}, ['^' currPrefix runs{sj, r}],kernel_size);
                            display([SJs{sj} ', ' runs{r} ' is done'])
                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                display(['Step 5, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                run_dir = fullfile(sesPath);
                                B5_smoothing(run_dir, SJs{sj}, ['^' currPrefix SJs{sj} '.*' runs{sj, r}(regexp(runs{sj,r}, 'run'):end)],kernel_size);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end
            currPrefix=['s' num2str(unique(kernel_size)) currPrefix];
        
        case 6 %% CompCorr
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                        counter = 0;
                        for r = 1:size(runs, 2)
                            display(['Step 6, CompCorr: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            data_dir = fullfile(funcPath, 'func');
                            if sj_space == 1
                                anat_dir = fullfile(funcPath, 'anat'); % iy_sub-004_ses-1_T1w

                                if counter == 0

                                    spm('defaults','fmri');
                                    spm_jobman('initcfg');

                                    warning off

                                    f1 = spm_select('List', anat_dir, '^iy.*\.nii$');
                                    numVols = size(f1,1);
                                    parameter_file = cellstr([repmat([anat_dir filesep], numVols, 1) f1 ]);

                                    f2 = spm_select('List', mni_wm_mask)
                                    f3 = spm_select('List', mni_csf_mask)
                                    numVols2 = size(f2,1);
                                    V=spm_vol([repmat([data_dir filesep], numVols2, 1) f2 ',:']);
                                    V2=spm_vol([repmat([data_dir filesep], numVols2, 1) f3 ',:']);

                                    for i=1:size(V,1)
                                        functional_imgs{i} = [repmat([data_dir filesep], numVols2, 1) f2 repmat([',' int2str(i)] , numVols2, 1)];
                                    end
                                    for i=(size(V,1)+1):size(V,1)*2
                                        functional_imgs{i} = [repmat([data_dir filesep], numVols2, 1) f3 repmat([',' int2str(i)] , numVols2, 1)];
                                    end

                                    matlabbatch{1}.spm.spatial.normalise.write.subj.def = parameter_file;
                                    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = functional_imgs';

                                    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;78 76 85];
                                    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                                    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

                                    spm_jobman('run', matlabbatch)
                                    clear jobs

                                    wm_mask=[anat_dir filesep 'sj_wm_mask_eroded.nii']; %white matter mask file
                                    csf_mask=[anat_dir filesep 'sj_csf_mask_eroded.nii'];

                                    movefile(['.../toolboxes' filesep 'wwm_mask_eroded.nii'],wm_mask);
                                    movefile(['.../toolboxes' filesep 'wcsf_mask_eroded.nii'],csf_mask);

                                    counter = 1;
                                end

                                wm_mask=[anat_dir filesep 'sj_wm_mask_eroded.nii']; %white matter mask file
                                csf_mask=[anat_dir filesep 'sj_csf_mask_eroded.nii'];

                            else
                                wm_mask = mni_wm_mask;
                                csf_mask = mni_csf_mask;
                            end
                            B6_compcorr(data_dir, SJs{sj}, ['^' currPrefix runs{sj, r}], numComp, wm_mask, csf_mask);
                        end
                    elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                        for ses = 1:sessNum
                            counter = 0;
                            for r = 1:size(runs, 2)
                                ses_dir = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                data_dir = fullfile(ses_dir, 'func');
                                if sj_space == 1
                                    anat_dir = fullfile(ses_dir, 'anat');

                                    if counter == 0

                                        spm('defaults','fmri');
                                        spm_jobman('initcfg');

                                        warning off

                                        delete(['.../toolboxes' filesep 'wwm_mask_eroded.nii']);
                                        delete(['.../toolboxes' filesep 'wcsf_mask_eroded.nii']);

                                        f1 = spm_select('List', anat_dir, '^iy.*\.nii$');
                                        parameter_file = cellstr([anat_dir filesep f1]);

                                        f2 = spm_select('List', '.../toolboxes', 'wm_mask_eroded.nii') 
                                        f3 = spm_select('List', '.../toolboxes', 'csf_mask_eroded.nii')
                                        functional_imgs = {['.../toolboxes' filesep f2], ['.../toolboxes' filesep f3]};

                                        matlabbatch{1}.spm.spatial.normalise.write.subj.def = parameter_file;
                                        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = functional_imgs';

                                        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;78 76 85];
                                        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                                        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

                                        f4 = spm_select('List', data_dir, '^rsub.*\.nii$');

                                        matlabbatch{2}.spm.spatial.realign.write.data = {
                                            [data_dir filesep f4(1,:) ',1']
                                            ['.../toolboxes' filesep 'wwm_mask_eroded.nii']
                                            ['.../toolboxes' filesep 'wcsf_mask_eroded.nii']
                                            };
                                        matlabbatch{2}.spm.spatial.realign.write.roptions.which = [1 0];
                                        matlabbatch{2}.spm.spatial.realign.write.roptions.interp = 4;
                                        matlabbatch{2}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
                                        matlabbatch{2}.spm.spatial.realign.write.roptions.mask = 1;
                                        matlabbatch{2}.spm.spatial.realign.write.roptions.prefix = 'r';

                                        spm_jobman('run', matlabbatch)
                                        clear jobs

                                        wm_mask=[anat_dir filesep 'sj_wm_mask_eroded.nii']; %white matter mask file
                                        csf_mask=[anat_dir filesep 'sj_csf_mask_eroded.nii'];

                                        delete(['.../toolboxes' filesep 'wwm_mask_eroded.nii']);
                                        delete(['.../toolboxes' filesep 'wcsf_mask_eroded.nii']);

                                        movefile(['.../toolboxes' filesep 'rwwm_mask_eroded.nii'],wm_mask);
                                        movefile(['.../toolboxes' filesep 'rwcsf_mask_eroded.nii'],csf_mask);

                                        counter = 1;
                                    end

                                    wm_mask=[anat_dir filesep 'sj_wm_mask_eroded.nii']; %white matter mask file
                                    csf_mask=[anat_dir filesep 'sj_csf_mask_eroded.nii'];
                                else
                                    wm_mask = mni_wm_mask;
                                    csf_mask = mni_csf_mask;
                                end
                                display(['Step 6, cc: ' SJs{sj} ', ' runs{sj, r}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                data_dir = fullfile(sesPath);
                                B6_compcorr(data_dir, SJs{sj}, ['^' currPrefix SJs{sj} '.*' runs{sj, r}(regexp(runs{sj,r}, 'run'):end)], numComp, wm_mask, csf_mask);
                            end
                        end
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                    end
                end
            end



        otherwise
            display('######################################################################################')
            display(['############################## Case ' num2str(n) ' does not exsist ##############################'])
            display('######################################################################################')
    end
    disp('Done.')

end