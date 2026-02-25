
%% INPUT ########################################################## 

% NOTE: before starting run 1, run the randomizer script or load
%	the assignment of cues into working memory
%       change the root_folder, and run vaiable
% remember to instruct participant accordingly!

subj = 25;	% subject number
sess = 2;	% session number
run = 7;	% run number

% dry run or real experiment?
stimulation_possible = 1; % 0 = no, 1 = yes
with_fmri = 1; % 0 = no fmri, 1 = with fmri

root_folder = ['...\IMACU\experiment']; % where pngs of cues are, etc

%%%#################################################################

%% specifics

% cues
preCue_time = 1; % in s
cue_time = 5; % e.g. 5s stimulation / imagery

% stim_cue determined in randomizer
% GET assignment from randomizer script beforehand

% colour assignment changes every run
run_colours = cue_colours(run, :);
run_colours_string = cue_colours_string(run, :);

% load according assignment cue
assignment_cue = fullfile(root_folder, 'cues', sprintf('leg_anterior_%d_%s_posteriorU_%d_%s.png', ...
point_label_assignment(1), run_colours_string{point_label_assignment(1)}, ...
point_label_assignment(2), run_colours_string{point_label_assignment(2)}));

% load this run colours check is stim cued by flipped or sqr

 if stim_cue == 1 % raute (flipped square)
        cue1 = fullfile(root_folder, 'cues', [run_colours_string{1} '_raute.png']); % stim1
        cue2 = fullfile(root_folder, 'cues', [run_colours_string{2} '_raute.png']); % stim2
        cue3 = fullfile(root_folder, 'cues', [run_colours_string{3} '_raute.png']); % stim3
        cue4 = fullfile(root_folder, 'cues', [run_colours_string{1} '_square.png']); % imag1
        cue5 = fullfile(root_folder, 'cues', [run_colours_string{2} '_square.png']); % imag2
        cue6 = fullfile(root_folder, 'cues', [run_colours_string{3} '_square.png']); % imag3

 elseif stim_cue == 0 % square
        cue1 = fullfile(root_folder, 'cues', [run_colours_string{1} '_square.png']); % stim1
        cue2 = fullfile(root_folder, 'cues', [run_colours_string{2} '_square.png']); % stim2
        cue3 = fullfile(root_folder, 'cues', [run_colours_string{3} '_square.png']); % stim3
        cue4 = fullfile(root_folder, 'cues', [run_colours_string{1} '_raute.png']); % imag1
        cue5 = fullfile(root_folder, 'cues', [run_colours_string{2} '_raute.png']); % imag2
        cue6 = fullfile(root_folder, 'cues', [run_colours_string{3} '_raute.png']); % imag3
 end

% read pngs
[assignment_cue, ~, ~] = imread(assignment_cue);

[cue1, ~, alpha_c1] = imread(cue1);
cue1(:, :, 4) = alpha_c1;
[cue2, ~, alpha_c2] = imread(cue2);
cue2(:, :, 4) = alpha_c2;
[cue3, ~, alpha_c3] = imread(cue3);
cue3(:, :, 4) = alpha_c3;
[cue4, ~, alpha_c4] = imread(cue4);
cue4(:, :, 4) = alpha_c4;
[cue5, ~, alpha_c5] = imread(cue5);
cue5(:, :, 4) = alpha_c5;
[cue6, ~, alpha_c6] = imread(cue6);
cue6(:, :, 4) = alpha_c6;

% load 4 possible versions of response cue (spiral)
[ans1, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver1.png'));
[ans2, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver2.png'));
[ans3, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver3.png'));
[ans4, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver4.png'));

%% stimulation
 
cond_names = {'stim point 1', 'stim point 2', 'stim point 3', 'imag point 1', 'imag point 2', 'imag point 3'};

if sess == 1 % all 6 conds

    condNum = 6;
    stim = [ones(5,1); ones(5,1)*2; ones(5,1)*3; ones(5,1)*4; ones(5,1)*5; ones(5,1)*6]; 

elseif sess == 2 % just imagery

    condNum = 3;
    stim = [ones(5,1)*4; ones(5,1)*5; ones(5,1)*6]; 

end

% randomize trial order
rando = randperm(length(stim));
stim_vec = [stim(rando)]; 
trialsum = length(stim_vec);

%% Psychtoolbox Setup
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Select the external screen if it is present, else revert to the native
% screen
screenNumber = max(screens);
if screenNumber >= 1
    screenNumber2 = screenNumber-1;
end

% supress tests and warnings
Screen('Preference', 'Verbosity', 0);
Screen('Preference', 'SkipSyncTests',1);
Screen('Preference', 'VisualDebugLevel',0);

% Open an on screen window
[window, windowRect] = Screen('OpenWindow', screenNumber, [0 0 0]);
% Set the text preferences
Screen('TextFont', window, 'Arial');
Screen('TextSize', window, 80); % adjust to situation
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
DrawFormattedText(window, 'waiting to start...', 'center', 'center', [255 255 255]);
Screen(window,'flip');

midImSq = [0 0 300 300]; % middle sized images, e.g. assignment cue
[midIm, xOffsetsigSb, yOffsetsigSb] = CenterRect(midImSq, windowRect);
tinyImSq = [0 0 50 50]; % very small images, e.g. pre-cues
[tinyIm, xOffsetsigSa, yOffsetsigSa] = CenterRect(tinyImSq, windowRect);
smImSq = [0 0 100 100]; % small images, e.g. cue
[smallIm, xOffsetsigS, yOffsetsigS] = CenterRect(smImSq, windowRect);

% make textures 
ascue_t = Screen('MakeTexture', window, assignment_cue);
cue1_t = Screen('MakeTexture', window, cue1);
cue2_t = Screen('MakeTexture', window, cue2);
cue3_t = Screen('MakeTexture', window, cue3);
cue4_t = Screen('MakeTexture', window, cue4);
cue5_t = Screen('MakeTexture', window, cue5);
cue6_t = Screen('MakeTexture', window, cue6);

ans1_t = Screen('MakeTexture', window, ans1);
ans2_t = Screen('MakeTexture', window, ans2);
ans3_t = Screen('MakeTexture', window, ans3);
ans4_t = Screen('MakeTexture', window, ans4);

%% prepare response options CHANGE

% draw up options: ? in the middle, sourrounded by the spiral of thickness 1-4,
% always in clockwise order but not always at same position: rotate answers
% so options are (gonna look awkward as a line of text):
%    ver1                  ver2                ver3                ver4
%     1                     4                   3                   2
%   4 ? 2                 3 ? 1               2 ? 4               1 ? 3
%     3                     2                   1                   4
% therefore buttons A, B, C and D should correspond to...
% A=1,B=4,C=2,D=3       A=4,B=3,C=1,D=2     A=3,B=2,C=4,D=1     A=2,B=1,C=3,D=4

% first as an array
button_order_options = [1, 4, 2, 3; ...
                        4, 3, 1, 2; ...
                        3, 2, 4, 1; ...
                        2, 1, 3, 4];
vers = 1:4;
% here also as an optional string to display
% ver1 = ' 1 \n \n 4        ?        2 \n \n 3 '; % for whatever reason, psychtoolbox doesn't accept \t as a tab
% ver2 = ' 4 \n \n 3        ?        1 \n \n 2 ';
% ver3 = ' 3 \n \n 2        ?        4 \n \n 1 ';
% ver4 = ' 2 \n \n 1        ?        3 \n \n 4 ';
% versions = {ver1, ver2, ver3, ver4}; % strings to pick
versions = {ans1_t, ans2_t, ans3_t, ans4_t}; % with images of spiral

% when to rotate?
rotation = zeros(1, trialsum);
rotation(1:round(trialsum/3)) = 1; 
rotation = rotation(randperm(trialsum));
% orientation per trial:
orientation = mod(cumsum(rotation), 4) + 1; % ~* modulo magic *~

%% Keyboard Tings

KbName('UnifyKeyNames')
escapeKey = KbName('ESCAPE');
%button box + alternatives in case I'm wrong (i was not)
buttonA = KbName('1');
buttonA2 = KbName('1!');
buttonB = KbName('2'); 
buttonB2 = KbName('2@');
buttonC = KbName('3'); 
buttonC2 = KbName('3#');
buttonD = KbName('4'); 
buttonD2 = KbName('4$');
%trigger
fmriTrig1 = KbName('5%'); 
fmriTrig2 = KbName('5');

% show assignment_cue 
Screen('DrawTexture', window, ascue_t, [], midIm);
Screen(window,'flip');

while 1
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown

        KbReleaseWait;

        %back to fixation cross after response
        DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
        Screen(window,'flip');

        break;
    end
end


%% ITIs

%#1: create a vector of 6 ITIs as many ITIs as trials, rep x conditions
ITI_vec = [2 3 3.5 4 5]; % 2 to 5 seconds ITIs
%#2: do permutations of this vector per condition (in case its just 2 conditions, theres double the trials so it should even out)
ITI_vec1 = [repmat(ITI_vec(randperm(length(ITI_vec))),[1,6])];
ITI_vec2 = [repmat(ITI_vec(randperm(length(ITI_vec))),[1,6])];
%#3: permate vector in order of trial types --> pre and post response ISIs
ITIs1 = [ITI_vec1(rando)]'; % three hard-coded ones up front
ITIs2 = [ITI_vec2(rando)]'; % three hard-coded ones up front


%% Design Matrix

% create design mat with place holders:
DesignMat = [   stim_vec               ...    %1. Trial Number
                zeros(trialsum,1)      ...    %2. Trial onset time(s)
                ITIs2                  ...    %3. pre ITI
                zeros(trialsum,1)      ...    %4. Response
                ITIs1 ];                      %5. post ITI for this trial

%% fMRI Trigger

t_start_mri = 0;
preTrial = 2; % time to wait from trigger to first trial

if with_fmri == 1
    activeKeys = [fmriTrig1, fmriTrig2];
    % RestrictKeysForKbCheck(activeKeys);
    clc;
    disp('Waiting for trigger...')

    while 1
        [keyIsDown,timeSecs,keyCode]=PsychHID('KbCheck');
    %     [keyIsDo5wn, timeSe5cs, keyCode] = KbCheck;
        keyCode = find(keyCode, 1);
        if keyIsDown        
            if keyCode == fmriTrig1 ||  keyCode == fmriTrig2
                t_start_mri = GetSecs;
                break;
            end
            % To condense multiple 'keyDown' events into a single event, 
            % we wait until all keys have been released.
            KbReleaseWait;
        end
    end
    disp('fMRI scanning started.')

else
    t_start_mri=GetSecs;
end

%% trial loop
% trials work as follows:
% 2s wait after trigger
% then start
% 3s before trial start cue to acupuncturist if applicable
% 1s pre-cue (small version of cue)
% 5s actual cue (and Stimulation if stim trial)
% 2-5s ITI: +
% 3s for button press 1 2 3 or 4 with options displayed
% 2-5s ITI again: +

%% sounds
% mapping to acupuncturist: position in point_label_assignment = her
% "number" for that point
% e.g. point_label_assignment = [3 1 2] -> YC hears one beep in
% condition 3, two beeps in condition 1 and three beeps in condition 2

yh1 = find(point_label_assignment == 1);
yh2 = find(point_label_assignment == 2);
yh3 = find(point_label_assignment == 3);

% audio cue takes about 0.25s so:
extra_y_time = 1.75; % time needed to prepare ahead of precue

Fs0 = 44100;
[y1, Fs1] = audioread(['cues\sound_cue' num2str(yh1) '.wav']);
[y2, Fs2] = audioread(['cues\sound_cue' num2str(yh2) '.wav']);
[y3, Fs3] = audioread(['cues\sound_cue' num2str(yh3) '.wav']);

InitializePsychSound;
panhandle = PsychPortAudio('Open', [], 1, 0, Fs0, 1);

%% actual trial loop
esc=0;
for t=1:trialsum

    if (t==1)
        disp('Waiting before initiating trial loop...')
        WaitSecs(preTrial); % 2 s 
    end

    condition = stim_vec(t);

    disp(['Trial ' num2str(t) ', Condition: ' cond_names{condition} ' (Colour: ' run_colours_string{mod(condition-1,3)+1} ')']);

    fprintf('ITI 2: %d\n', ITIs2(t));
    WaitSecs(ITIs2(t)-extra_y_time);
    
    switch condition
        case 1 % stim 1

            PsychPortAudio('FillBuffer', panhandle, y1');
            PsychPortAudio('Start', panhandle, 1, 0, 1);

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            Screen('DrawTexture', window, cue1_t, [], tinyIm);
            Screen(window,'flip');

            WaitSecs(preCue_time);
            
            % cue1
            Screen('DrawTexture', window, cue1_t, [], smallIm);
            Screen(window,'flip');
                 
            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
        case 2 % stim 2 
            % precue

            PsychPortAudio('FillBuffer', panhandle, y2');
            PsychPortAudio('Start', panhandle, 1, 0, 1);

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            Screen('DrawTexture', window, cue2_t, [], tinyIm);
            Screen(window,'flip');

            WaitSecs(preCue_time);
            
            % cue2
            Screen('DrawTexture', window, cue2_t, [], smallIm);
            Screen(window,'flip');
                        
            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
        case 3 % stim 3 

            PsychPortAudio('FillBuffer', panhandle, y3');
            PsychPortAudio('Start', panhandle, 1, 0, 1);

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            % precue
            Screen('DrawTexture', window, cue3_t, [], tinyIm);
            Screen(window,'flip');

            WaitSecs(preCue_time);
            
            % cue3
            Screen('DrawTexture', window, cue3_t, [], smallIm);
            Screen(window,'flip');

            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
        case 4 % imag 1

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            % precue
            Screen('DrawTexture', window, cue4_t, [], tinyIm);
            Screen(window,'flip');
            
            WaitSecs(preCue_time);
            
            % cue4
            Screen('DrawTexture', window, cue4_t, [], smallIm);
            Screen(window,'flip');
            
            WaitSecs(cue_time);
            
        case 5 % imag 2 

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            % precue
            Screen('DrawTexture', window, cue5_t, [], tinyIm);
            Screen(window,'flip');
            
            WaitSecs(preCue_time);
            
            % cue5 no stim
            Screen('DrawTexture', window, cue5_t, [], smallIm);
            Screen(window,'flip');

            WaitSecs(cue_time);


        case 6 % imag 3 

            WaitSecs(extra_y_time);
            trial_start = GetSecs-t_start_mri;

            % precue
            Screen('DrawTexture', window, cue6_t, [], tinyIm);            
            Screen(window,'flip');

            WaitSecs(preCue_time);

            % cue6 no stim
            Screen('DrawTexture', window, cue6_t, [], smallIm);            
            Screen(window,'flip');

            WaitSecs(cue_time);

    end

    % back to fixation cross
    DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
    Screen(window,'flip');
    
    %Wait duration of the current ITI before response
    fprintf('ITI 1: %d\n', ITIs1(t));
    WaitSecs(ITIs1(t));

    % RESPONSE
    fprintf(['version: ', num2str(orientation(t))]);
    % DrawFormattedText(window, versions{orientation(t)}, 'center', 'center', [255 255 255]);
    Screen('DrawTexture', window, versions{orientation(t)}, [], midIm);
    Screen(window,'flip');

    resp_start = GetSecs;
    tnow = GetSecs;
    pressed = 0;
    
    % don't let the MRI answer
    RestrictKeysForKbCheck([escapeKey, buttonA, buttonA2, buttonB, buttonB2, buttonC, buttonC2, buttonD, buttonD2]);
    
    while tnow < resp_start+3 % 3s for answer 
        tnow = GetSecs;
        [keyIsDown, ~, keyCode] = KbCheck;

        if keyCode(KbName('ESCAPE')) == 1
            esc=1;
            break;
        end

        keyCode = find(keyCode, 1);
        if keyIsDown        
            if keyCode == buttonA || keyCode == buttonA2
                pressed = 3;
            elseif keyCode == buttonB || keyCode == buttonB2
                pressed = 4;
            elseif keyCode == buttonC || keyCode == buttonC2
                pressed = 2;
            elseif keyCode == buttonD || keyCode == buttonD2
                pressed = 1;
            end
            KbReleaseWait;
            
            %back to fixation cross after response
            DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
            Screen(window,'flip');
        end
    end
  
    % record response
    if pressed > 0
        response = button_order_options(orientation(t), pressed);
    else
        response = 0;
    end

    %back to fixation cross anyhow
    DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
    Screen(window,'flip');
       
    fprintf(['response: ', num2str(response), '\n']);

    if esc==1 % experimentor panic button
        disp('execution stopped by user --> ESC');
        sca
        break;
    end

    % fill design mat
    DesignMat(t, 2) = trial_start;
    DesignMat(t, 4) = response;
end

PsychPortAudio('Close', panhandle);

clear screen

%% write log 

%Get datetime
dt=datestr(now,'dd.mm.YYYY_HH-MM-SS');
save(['Logs\DesignMat_subject_' num2str(subj) '_sess_' num2str(sess) '_run_' num2str(run) '_' dt '.mat'],"DesignMat")
DesignMat_C = arrayfun(@num2str,DesignMat,'UniformOutput',false);

% Creating the log file to be filled during the experiment
fileID = fopen(['Logs\LogFile_subject_' num2str(subj) '_sess_' num2str(sess)  '_run_' num2str(run) '_' dt '.tsv'],'w');       %creates the empty logfile in the directory you specify 
formatSpec = '%s\t%s\t%s\t%s\t%s\n';        %specifies the format of each variable that you put into the logfile
labels = ["TrialNr" "onset" "preITI" "Response" "postITI"];     %Specifies labels for each column that should be in the logfile.
fprintf(fileID, formatSpec, labels);        %This actially prints the line labels into the first row of the log file.

%%% Writing log file in .tsv format
for i=1:trialsum
    log = DesignMat_C{i,1:5};     %define the next row of the log file. This is done after the end of each trial to save the actial values. 
    fprintf(fileID, formatSpec, DesignMat_C{i,1:5});       %print the previously defined row in the log file.
end
% Closing the log file
fclose(fileID);     % This closes the log file at the end of your experiment.

disp('Done! :)');