%% decoding BATCH %%

% required toolboxes:
% the decoding toolbox TDT
%  https://doi.org/10.3389/fninf.2014.00088

%#####################################################
%#################### INPUT ##########################
%#####################################################

clear runs


%SPM-path
SPM_path  = '.../toolboxes/spm12';

%data source directory
src_dir      = '.../IMACU/Data';

addpath(genpath('.../IMACU/C_statistical_analysis-main'));
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
display('Subjects to exclude:')
excludeSJ = [11] % exclude those sjs: 1 22 | 9 16 25 30 | 8 10 23 | 4 17 18 24

%% selection of analysis steps (1-5) to be performed
analysis_switch = [ 1 2 ]; % 1 2 3 4 5

zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files)
    for z = 1:size(zip_files, 1)
        gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1:  2nd level flexfact on univar data
con_images=1:3; 
session_2nd = 1;
outputfolder_2nd_b = '2nd_level_N25_sess2_flexfact'; % folder that will contain the created SPM file
dir_1st_b   = {'ses-2/1st_level_C_s8wr_hmcc'}; % Name of corresponding first Level Analysis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  2nd level flexible factorial over confusion results
pref2 = 's5w';
conds2 = {'TrueAcu', 'AcuToC1', 'AcuToC2', 'C1ToAcu', 'TrueCp1', 'Cp1ToC2', 'C2ToAcu', 'Cp2ToC1', 'TrueCp2'};
outputfolder_2nd2 = {'2nd_level_DEC_N25_sess1_to_sess2_imag_FlexFact_confusion'}; % folder that will contain the created SPM file
folder2_inst1 = ['DEC']; % 'ses-1' filesep if in first session
dir_1st2   = {'D_ImagAcu_ImagC1_ImagC2_ses_1_to_2'}; % Name of corresponding first Level Analysis %  e.g. D4_StimAcu_StimC1_StimC2 D4_ImagAcu_ImagC1_ImagC2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        %% 2nd level flexible factorial univar
        case 1

            SJin=SJs;
            SJin(excludeSJ)=[];
            E1_2ndLevel_FlexFact(src_dir, SJin, outputfolder_2nd_b, dir_1st_b, con_images, session_2nd);


        %% 2nd level flexible factorial multivar on confusion
        case 8

            SJin=SJs;
            SJin(excludeSJ)=[];
            for folds = 1:length(outputfolder_2nd2)
                this_2nd = outputfolder_2nd2{folds};
                this_1st = dir_1st2{folds};
                E2_2ndLevel_FlexFact_confusion(src_dir, SJin, this_2nd, [folder2_inst1 filesep this_1st], conds2, pref2);
            end

    end
end
disp('Done.')

