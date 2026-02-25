% test run

% very similar to E experimental script --> datailed documentation there

sess = 1;
run = 1;

root_folder = ['...\IMACU\experiment\'];


%% specifics

% cues
preCue_time = 1; %in s
cue_time = 5; % 1 for testing, 5 for real

for r = 1
    this_rand_order = randperm(length(add_cue));
    rand_order = [rand_order; this_rand_order];
    cue_colours = [cue_colours; add_cue(this_rand_order)];
    cue_colours_string = [cue_colours_string; add_strings(this_rand_order)];
    disp(['Colours for run ' num2str(r) ': ' add_strings{this_rand_order(1)} ', ' add_strings{this_rand_order(2)} ', ' add_strings{this_rand_order(3)}])
end

%% colour assignment
run_colours = cue_colours(run, :);
run_colours_string = cue_colours_string(run, :);

%% assignment cue
assignment_cue = fullfile(root_folder, 'cues', sprintf('leg_anterior_%d_%s_posteriorU_%d_%s.png', ...
    point_label_assignment(1), run_colours_string{point_label_assignment(1)}, ...
    point_label_assignment(2), run_colours_string{point_label_assignment(2)}));

% stim cued by flipped or upright square

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

% load pngs witz transparent background
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

% answer images
[ans1, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver1.png'));
[ans2, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver2.png'));
[ans3, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver3.png'));
[ans4, ~, ~] = imread(fullfile(root_folder, 'cues', 'ver4.png'));

%% stimulation rand and order

cond_names = {'stim point 1', 'stim point 2', 'stim point 3', 'imag point 1', 'imag point 2', 'imag point 3'};

if sess == 1 % all 6 conds,
    % changed to 2 so we have 2 rand trials per cond
    %and changed to stim1 so stim is the final order (points123 then rand )
    condNum = 6;
    stim1 = [ones(2,1); ones(2,1)*2; ones(2,1)*3; ones(2,1)*4; ones(2,1)*5; ones(2,1)*6];

elseif sess == 2 % jus ?

    condNum = 3;
    stim1 = [ones(10,1)*4; ones(10,1)*5; ones(10,1)*6];

end

% check later - trialsum is needed in loop; changed for test
rando = randperm(length(stim1));
stim_vec1 = stim1(rando);

trials1 = [1; 4; 2; 5; 3; 6; 1; 4; 2; 5; 3; 6];
stim_vec2 = [trials1; stim_vec1];

trialsum = length(stim_vec2);
rando = randperm(trialsum);

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

% draw buttons for visual feedback while answering
button_radius = 25;
button_upper =  [round(windowRect(3)/2)-button_radius   yOffsetsigSb-button_radius*3            round(windowRect(3)/2)+button_radius        yOffsetsigSb-button_radius];
button_left =   [xOffsetsigSb-button_radius*3           round(windowRect(4)/2)-button_radius    xOffsetsigSb-button_radius                  round(windowRect(4)/2)+button_radius];
button_lower =  [round(windowRect(3)/2)-button_radius   yOffsetsigSb+midImSq(4)+button_radius   round(windowRect(3)/2)+button_radius        yOffsetsigSb+midImSq(4)+button_radius*3];
button_right =  [xOffsetsigSb+midImSq(3)+button_radius  round(windowRect(4)/2)-button_radius    xOffsetsigSb+midImSq(3)+button_radius*3     round(windowRect(4)/2)+button_radius];
all_buttons = [button_upper' button_left' button_lower' button_right'];

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

%% response options

button_order_options = [1, 4, 2, 3; ...
    4, 3, 1, 2; ...
    3, 2, 4, 1; ...
    2, 1, 3, 4];
vers = 1:4;
% here also as a string to display
% ver1 = ' 1 \n \n 4        ?        2 \n \n 3 '; % for whatever reason, psychtoolbox doesn't accept \t as a tab
% ver2 = ' 4 \n \n 3        ?        1 \n \n 2 ';
% ver3 = ' 3 \n \n 2        ?        4 \n \n 1 ';
% ver4 = ' 2 \n \n 1        ?        3 \n \n 4 ';
% versions = {ver1, ver2, ver3, ver4}; % strings to pick
versions = {ans1_t, ans2_t, ans3_t, ans4_t}; % pngs of spiral

% when to rotate?
rotation = zeros(1, trialsum);
rotation(1:round(trialsum/3)) = 1; 
rotation = rotation(randperm(trialsum));
% orientation per trial:
orientation = mod(cumsum(rotation), 4) + 1; % ~* modulo magic *~


%% Keyboard Tings

KbName('UnifyKeyNames')
escapeKey = KbName('ESCAPE');
%button box + alternatives in case I'm wrong (I'm not)
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
% change added _2 so diff names then actual exp runs
%#1: create a vector of 6 ITIs as many ITIs as trials, rep x conditions
ITI_vec_2 = [2 2.5 3 3.5 4 5]; % 2 to 5 seconds ITIs CHANGED bc 18 trials sum
%#2: do permutations of this vector per condition (in case its just 2 conditions, theres double the trials so it should even out)
ITI_vec1_1 = [repmat(ITI_vec_2(randperm(length(ITI_vec_2))),[1,4])];
ITI_vec2_2 = [repmat(ITI_vec_2(randperm(length(ITI_vec_2))),[1,4])];
%#3: permate vector in order of trial types --> pre and post response ISIs
ITIs1 = [ITI_vec1_1(rando)]';
ITIs2 = [ITI_vec2_2(rando)]';


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


%% actual trial loop
esc=0;
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

            % precue: small version of cue
            Screen('DrawTexture', window, cue1_t, [], tinyIm);
            Screen(window,'flip');
            
            WaitSecs(preCue_time);

            % cue1: hand with needle
            Screen('DrawTexture', window, cue1_t, [], smallIm);
            Screen(window,'flip');

            WaitSecs(cue_time);
            PsychPortAudio('Stop', panhandle, 1);

        case 2 % stim 2 

            PsychPortAudio('FillBuffer', panhandle, y2');
            PsychPortAudio('Start', panhandle, 1, 0, 1);

            WaitSecs(extra_y_time);

            % precue
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
            Screen('DrawTexture', window, cue4_t, [], tinyIm);
            Screen(window,'flip');

            WaitSecs(preCue_time);

            % cue4
            Screen('DrawTexture', window, cue4_t, [], smallIm);
            Screen(window,'flip');

            WaitSecs(cue_time);

        case 5 % imag 2 
            WaitSecs(extra_y_time);
            Screen('DrawTexture', window, cue5_t, [], tinyIm);
            Screen(window,'flip');

            WaitSecs(preCue_time);

            % cue5 no stim
            Screen('DrawTexture', window, cue5_t, [], smallIm);
            Screen(window,'flip');

            WaitSecs(cue_time);


        case 6 % imag 3 
            WaitSecs(extra_y_time);
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

    % RESPONSE + visual feedback of button press
    fprintf(['version: ', num2str(orientation(t))]);
    Screen('DrawTexture', window, versions{orientation(t)}, [], midIm);
    Screen('FillOval', window, [150], all_buttons)
    Screen(window,'flip');

    resp_start = GetSecs;
    tnow = GetSecs;
    pressed = 0;

    RestrictKeysForKbCheck([escapeKey, buttonA, buttonA2, buttonB, buttonB2, buttonC, buttonC2, buttonD, buttonD2]);

    while tnow < resp_start+3
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

                Screen('FillOval', window, [150 150 150 50; 150 150 150 50; 150 150 150 50], all_buttons)

            elseif keyCode == buttonB || keyCode == buttonB2

                pressed = 4;

                Screen('FillOval', window, [150 150 50 150; 150 150 50 150; 150 150 50 150], all_buttons)
            
            elseif keyCode == buttonC || keyCode == buttonC2

                pressed = 2;

                Screen('FillOval', window, [150 50 150 150; 150 50 150 150; 150 50 150 150], all_buttons)
            
            elseif keyCode == buttonD || keyCode == buttonD2

                pressed = 1;

                Screen('FillOval', window, [50 150 150 150; 50 150 150 150; 50 150 150 150], all_buttons)
            
            end
            Screen(window,'flip');
            KbReleaseWait;
            
            WaitSecs(1);

            %back to fixation cross after response
            
            DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
            Screen(window,'flip');
        end
    end

    if pressed > 0
        response = button_order_options(orientation(t), pressed);
    else
        response = 0;
    end

    %back to fixation cross anyhow
    DrawFormattedText(window, '+', 'center', 'center', [255 255 255]);
    Screen(window,'flip');

    fprintf(['response: ', num2str(response), '\n']);

    if esc==1
        disp('execution stopped by user --> ESC');
        sca
        break;
    end

end

WaitSecs(2);

PsychPortAudio('Close', panhandle);

clear screen









