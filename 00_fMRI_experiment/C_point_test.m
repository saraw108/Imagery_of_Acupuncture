% run or load randomizer before
% similar to E experimental script --> detailed documentation there

sess = 1;
run = 1;

root_folder = ['...\IMACU\experiment\'];

preCue_time = 1; 
cue_time = 5;

for r = 1
    this_rand_order = randperm(length(add_cue));
    rand_order = [rand_order; this_rand_order];
    cue_colours = [cue_colours; add_cue(this_rand_order)];
    cue_colours_string = [cue_colours_string; add_strings(this_rand_order)];
    disp(['Colours for run ' num2str(r) ': ' add_strings{this_rand_order(1)} ', ' add_strings{this_rand_order(2)} ', ' add_strings{this_rand_order(3)}])
end

if mod(sj_id-1, 4) < 2
    disp('Stimulation is cued by the flipped square!');
    stim_cue = 1;
elseif mod(sj_id-1, 4) > 1
    disp('Stimulation is cued by the upright square!');
    stim_cue = 0;
end


run_colours = cue_colours(run, :);
run_colours_string = cue_colours_string(run, :);

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

%% stimulation rand and order

cond_names = {'stim point 1', 'stim point 2', 'stim point 3'}; %changed
condNum = 3;
stim1 = [ones(2,1); ones(2,1)*2; ones(2,1)*3];

stim_vec1 = stim1;
stim_vec2 = [1;2;3;1;2;3]; %used in loop

trialsum = length(stim_vec2);



%% Psychtoolbox
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
% Set the text stuff
Screen('TextFont', window, 'Arial');
Screen('TextSize', window, 80); % adjust if it looks stupid
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
DrawFormattedText(window, 'waiting to start...', 'center', 'center', [255 255 255]);
Screen(window,'flip');

midImSq = [0 0 300 300];
[midIm, xOffsetsigSb, yOffsetsigSb] = CenterRect(midImSq, windowRect);
tinyImSq = [0 0 50 50];
[tinyIm, xOffsetsigSa, yOffsetsigSa] = CenterRect(tinyImSq, windowRect);
smImSq = [0 0 100 100];
[smallIm, xOffsetsigS, yOffsetsigS] = CenterRect(smImSq, windowRect);

button_radius = 25;
button_upper =  [round(windowRect(3)/2)-button_radius   yOffsetsigSb-button_radius*3            round(windowRect(3)/2)+button_radius        yOffsetsigSb-button_radius];
button_left =   [xOffsetsigSb-button_radius*3           round(windowRect(4)/2)-button_radius    xOffsetsigSb-button_radius                  round(windowRect(4)/2)+button_radius];
button_lower =  [round(windowRect(3)/2)-button_radius   yOffsetsigSb+midImSq(4)+button_radius   round(windowRect(3)/2)+button_radius        yOffsetsigSb+midImSq(4)+button_radius*3];
button_right =  [xOffsetsigSb+midImSq(3)+button_radius  round(windowRect(4)/2)-button_radius    xOffsetsigSb+midImSq(3)+button_radius*3     round(windowRect(4)/2)+button_radius];
all_buttons = [button_upper' button_left' button_lower' button_right'];

% make textures
%ascue_t = Screen('MakeTexture', window, assignment_cue);



%% Keyboard Tings

KbName('UnifyKeyNames')
escapeKey = KbName('ESCAPE');
%button box + alternatives in case I'm wrong
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
%Screen('DrawTexture', window, ascue_t, [], midIm);
%Screen(window,'flip');

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

ITIs1 = [3 3 3 3 3 3]; % hardcode for test run
ITIs2 = [3 3 3 3 3 3];


%% sounds

% mapping to acupuncturist: position in point_label_assignment = her
% "number" for that point
% e.g. point_label_assignment = [3 1 2] -> Yinghui hears one beep in
% condition 3, two beeps in condition 1 and three beeps in condition 2

yh1 = find(point_label_assignment == 1);
yh2 = find(point_label_assignment == 2);
yh3 = find(point_label_assignment == 3);

extra_y_time = 1.75;

Fs1 = 44100;
[y1, Fs1] = audioread(['cues\sound_cue' num2str(yh1) '.wav']);
[y2, Fs2] = audioread(['cues\sound_cue' num2str(yh2) '.wav']);
[y3, Fs3] = audioread(['cues\sound_cue' num2str(yh3) '.wav']);
InitializePsychSound;
panhandle = PsychPortAudio('Open', [], 1, 0, Fs1, 1);

%% part 1: stim every point x2 so participants get familiar

esc=0;
for t=1:trialsum
    
    condition = stim_vec2(t);
    
    disp(['Trial ' num2str(t) ', Condition: ' cond_names{condition} ' (Colour: ' run_colours_string{mod(condition-1,3)+1} ')']);
    
    fprintf('ITI 2: %d\n', ITIs2(t));
    WaitSecs(ITIs2(t)-extra_y_time);
    
    switch condition
        case 1 % stim 1
            
            PsychPortAudio('FillBuffer', panhandle, y1');
            PsychPortAudio('Start', panhandle, 1, 0, 1);
            
            WaitSecs(extra_y_time);
            
            % cue1: hand with needle
            DrawFormattedText(window, 'This is point 1', 'center', 'center', [255 255 255]);
            Screen(window,'flip'); 
            
            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
        case 2 % stim 2 

            PsychPortAudio('FillBuffer', panhandle, y2');
            PsychPortAudio('Start', panhandle, 1, 0, 1);
            
            WaitSecs(extra_y_time);
            
            % cue2
            DrawFormattedText(window, 'This is point 2', 'center', 'center', [255 255 255]);
            Screen(window,'flip');
            
            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
        case 3 % stim 3 
            
            PsychPortAudio('FillBuffer', panhandle, y3');
            PsychPortAudio('Start', panhandle, 1, 0, 1);
            
            WaitSecs(extra_y_time);
            
            % cue3
            DrawFormattedText(window, 'This is point 3', 'center', 'center', [255 255 255]);
            Screen(window,'flip');
            
            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);
            
    end
    
    % back to fixation cross
    DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
    Screen(window,'flip');
    
    %Wait duration of the current ITI before response
    fprintf('ITI 1: %d\n', ITIs1(t));
    WaitSecs(ITIs1(t));
    
    resp_start = GetSecs;
    tnow = GetSecs;
    pressed = 0;
    
    %back to fixation cross anyhow
    DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
    Screen(window,'flip');
    
    if esc==1
        disp('execution stopped by user --> ESC');
        sca
        break;
    end
    
end


WaitSecs(2);

%PsychPortAudio('Close', panhandle);

%% part 2 - point recognition test for participant

rando = randperm(length(stim1));

stim_vec1 = stim1(rando);
stim_vec2 = stim_vec1;

trialsum = length(stim_vec2);

esc=0;

%loop for potential repetition

repy = 1;
cont = true;

while cont
    
    for t=1:trialsum
        
        condition = stim_vec2(t);
        
        disp(['Trial ' num2str(t) ', Condition: ' cond_names{condition} ' (Colour: ' run_colours_string{mod(condition-1,3)+1} ')']);
        
        % second ITI after response
        fprintf('ITI 2: %d\n', ITIs2(t));
        WaitSecs(ITIs2(t)-extra_y_time);
        
        switch condition
            case 1 % stim 1
                
                PsychPortAudio('FillBuffer', panhandle, y1');
                PsychPortAudio('Start', panhandle, 1, 0, 1);
                
                WaitSecs(extra_y_time);
                
                % cue1: hand with needle
                DrawFormattedText(window, 'Which point is that?', 'center', 'center', [255 255 255]);
                Screen(window,'flip'); %needed?
                
                WaitSecs(cue_time);
                PsychPortAudio('Stop', panhandle, 1);
                
            case 2 % stim 2 
                
                PsychPortAudio('FillBuffer', panhandle, y2');
                PsychPortAudio('Start', panhandle, 1, 0, 1);
                
                WaitSecs(extra_y_time);
                
                % cue2
                DrawFormattedText(window, 'Which point is that?', 'center', 'center', [255 255 255]);
                Screen(window,'flip');
                
                WaitSecs(cue_time);
                PsychPortAudio('Stop', panhandle, 1);
                
            case 3 % stim 3 
                
                PsychPortAudio('FillBuffer', panhandle, y3');
                PsychPortAudio('Start', panhandle, 1, 0, 1);
                
                WaitSecs(extra_y_time);
                
                % cue3
                DrawFormattedText(window, 'Which point is that?', 'center', 'center', [255 255 255]);
                Screen(window,'flip');
                
                WaitSecs(cue_time);
                PsychPortAudio('Stop', panhandle, 1);
                
        end
        
        % back to fixation cross
        DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
        Screen(window,'flip');
        
        %Wait duration of the current ITI before response
        fprintf('ITI 1: %d\n', ITIs1(t));
        WaitSecs(ITIs1(t));
        
        resp_start = GetSecs;
        tnow = GetSecs;
        pressed = 0;
        
        %back to fixation cross anyhow
        DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
        Screen(window,'flip');
        
        if esc==1
            disp('execution stopped by user --> ESC');
            sca
            break;
        end
        
        w = 0;
        while w ~= 1
            w = input('Did they answer? (1 = Yes, 0 = No): ');
        end
        
    end
    
    rep_yn = input('Do you want to repeat? (1 = Yes, 0 = No): ');
    
    if rep_yn == 1
        cont == true; %
    else
        cont == false; % exit loop
        break
    end
end

WaitSecs(2);

PsychPortAudio('Close', panhandle);

