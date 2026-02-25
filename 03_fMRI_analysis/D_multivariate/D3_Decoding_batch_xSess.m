function D3_Decoding_batch_xSess(beta_dirs, output_dir, labelnames_train, labelnames_test, from_to, dec_type, rad, masks, map_type, grid_C)
% This Batch Script first specifies what features should be decoded and
% executes the Decoding Cross Classification.
% The resulting accuracy maps should later be normalized and smoothed

% CREATE THE OUTPUT DIRECTORY IF IT DOES NOT EXIST YET
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% CONFIGURE THE DECODING TOOLBOX

clear cfg
cfg = decoding_defaults;
cfg.software = 'SPM12';

% Specify where the results should be saved
cfg.results.overwrite     = 1;
cfg.results.dir           = output_dir;

% DATA SCALING
cfg.scale.method          = 'z'; %z-score for all voxels over all samples (2 conditions x 6 runs =12)
cfg.scale.estimation      = 'all'; % only training, only test data, or all

% SEARCHLIGHT SPECIFICATIONS
cfg.analysis              = dec_type;   % or 'roi'
cfg.searchlight.unit      = 'voxels'; % or mm 
cfg.searchlight.radius    = rad; % 4 voxel or mm
cfg.searchlight.spherical = 0;  % only useful for mm
% The amount of information you want to have printed on the screen
% 0: no, 1: normal, output, 2: all)
cfg.verbose               = 0;  

% Method and model parameters 
cfg.decoding.method = 'classification';
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% OUTPUTS SPECIFICATION
cfg.results.output = map_type;
% 'accuracy', 'accuracy_minus_chance', 'sensitivity_minus_chance', 'specificity_minus_chance', 'balanced_accuracy_minus_chance', 'confusion_matrix', 'AUC_minus_chance'

% DISPLAY:
cfg.plot_selected_voxels  = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% This is by default set to 1, but if you repeat the same design again and again, it can get annoying...
cfg.plot_design           = 0;

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")

% cfg.files.mask = fullfile(beta_dir, 'mask.nii');

if dec_type == 'roi'
    cfg.files.mask = masks;
else
    cfg.files.mask = [beta_dirs(1,:) filesep 'mask.nii']; %lets make this simple
end

if grid_C == 1
 	cfg.parameter_selection.method = 'grid';
 	cfg.parameter_selection.parameter_range = {[0.0001 0.001 0.01 0.1 1 10 100 1000]};
 	cfg.parameter_selection.design.function.name = 'make_design_cv';
end

%% sort data
beta_dir_train = beta_dirs(from_to(1), :);
beta_dir_test = beta_dirs(from_to(2), :);

regressor_names_train = design_from_spm(beta_dir_train);
regressor_names_test = design_from_spm(beta_dir_test);

train1=find(strcmp(regressor_names_train(1,:),labelnames_train(1,1)));
train2=find(strcmp(regressor_names_train(1,:),labelnames_train(1,2)));
train3=find(strcmp(regressor_names_train(1,:),labelnames_train(1,3)));

test1=find(strcmp(regressor_names_test(1,:),labelnames_test(1,1)));
test2=find(strcmp(regressor_names_test(1,:),labelnames_test(1,2)));
test3=find(strcmp(regressor_names_test(1,:),labelnames_test(1,3)));

%Select Files by using a filter
f= {};
filename='beta_0000';
c       = []; % chunk
l       = []; % label
xc      = []; % xclass

for i = 1:length(train1)
    temp=num2str(train1(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_train), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 1;
    xc(length(xc)+1) = 1;
    
    temp=num2str(train2(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_train), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 2;
    xc(length(xc)+1) = 1;

    temp=num2str(train3(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_train), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 3;
    xc(length(xc)+1) = 1;

end

for i = 1:length(test1)
    temp=num2str(test1(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_test), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 1;
    xc(length(xc)+1) = 2;
    
    temp=num2str(test2(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_test), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 2;
    xc(length(xc)+1) = 2;

    temp=num2str(test3(i));
    f{length(f)+1} = spm_select('FPList',fullfile(beta_dir_test), [filename(1:end-length(temp)) temp '.(nii|img)']);
    c(length(c)+1)  = i;
    l(length(l)+1)  = 3;
    xc(length(xc)+1) = 2;

end

cfg.files.chunk  = c';
cfg.files.label  = l';
cfg.files.xclass = xc'; % 1: training, 2: test data
cfg.files.twoway = 0;   % 1 for average of  both tests (xclass "1 vs 2" and "2 vs 1")
cfg.files.name   = f';

% === Automatic Creation ===
% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_xclass_cv(cfg);

%$ Run DECODING
results = decoding(cfg);

