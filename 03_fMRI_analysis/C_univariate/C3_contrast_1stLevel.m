function C3_contrast_1stLevel(subj_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs, is_2nd, c_type)

%% Allocate SPM.mat file
if is_2nd == 1
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(analysisfolder, 'SPM.mat')};
else
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(subj_dir, analysisfolder, 'SPM.mat')};
end
    
% Number of Runs
% rundir=dir([subj_dir filesep 'func']);
numRuns=size(runs,2);       
% Number of contrasts
numCons  = numel(cnames);

% Cycle over contrast specifications
for c = 1:numCons
    if is_2nd == 1
        convec=cvecs(c, :);
    else
        convec=cvecs{c};
    end
    if n_hm>0
        convec=[convec zeros(1,n_hm)];
    end
    if n_cc>0
        convec=[convec zeros(1,n_cc)];
    end

    if strcmp(c_type, 'F')

        % Allocate t-contrast structure
        matlabbatch{1}.spm.stats.con.consess{c}.fcon.name    = cnames{c};
        if is_2nd == 1
            matlabbatch{1}.spm.stats.con.consess{c}.fcon.weights  = convec;
        else
            matlabbatch{1}.spm.stats.con.consess{c}.fcon.weights  = repmat(convec, 1, numRuns);
        end
        matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';

    elseif strcmp(c_type, 't')

        % Allocate t-contrast structure
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name    = cnames{c};
        if is_2nd == 1
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights  = convec;
        else
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights  = repmat(convec, 1, numRuns);
        end
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';

    end
end 

% Delete existing contrasts (1=yes)
matlabbatch{1}.spm.stats.con.delete = del;
    
% Run the job
fprintf(['Computing 1st Level Contrasts\n'])
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch