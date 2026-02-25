% run once, save sound files for ever
% writing sound files for IMACU

sf = 44100; % Sampling frequency
tPc = 0.05; % duration of cue
tCue = 5; % stim duration

mediumf = 300; % low freq 300 Hz

waveC = 0:1/sf:tPc;
waveC2 = 0:1/sf:tPc*2;

mediumSoundpc = sin(2 * pi * mediumf * waveC); % Cue 1 sound
noSound = zeros(1, tPc*4*sf); % filler gap
noSound2 = zeros(1, tPc*sf); % small gap
noSound3 = zeros(1, tPc*35*sf); % end gap

waveProper = 0:1/sf:tCue; 

mediumProper = sin(2 * pi * mediumf * waveProper); % Cue 1 sound

tPcf = linspace(0, tPc, sf * tPc);

sound_cue1 = [mediumSoundpc repelem(noSound2,3) repelem(noSound,4) noSound3 mediumProper]; 					% one beep
sound_cue2 = [mediumSoundpc noSound mediumSoundpc repelem(noSound2,2) repelem(noSound,3) noSound3 mediumProper];		% two beeps
sound_cue3 = [mediumSoundpc noSound mediumSoundpc noSound mediumSoundpc noSound2 repelem(noSound,2) noSound3 mediumProper];	% three beeps

audiowrite('Cues\sound_cue1.wav', sound_cue1, 44100);
audiowrite('Cues\sound_cue2.wav', sound_cue2, 44100);
audiowrite('Cues\sound_cue3.wav', sound_cue3, 44100);




