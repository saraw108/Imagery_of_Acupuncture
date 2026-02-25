% randomizer script

%% sj-ID
sj_id = 25;

point_lables = [1 2 3]; % position left: anterior (acu), middle: posterior upper, right: posterior lower 
add_cue = [2 1 0]; % colour of point 1, 2 and 3, irrespective of position on leg!!!
add_strings = {'yellow', 'blue', 'green'}; % colour of point 1, 2 and 3, irrespective of position on leg!!!

% randomize point numbering
for i = 1:sj_id
    point_label_assignment = point_lables(randperm(length(point_lables)));
end
disp(['Points for Subject ' num2str(sj_id) ': anterior point has label ' num2str(point_label_assignment(1)) ', posterior upper point has label ' num2str(point_label_assignment(2)) ', posterior lower point has label ' num2str(point_label_assignment(3)) ])

% randomize point colours per run
cue_colours_string = {};
cue_colours = [];
rand_order = [];
for r = 1:7
    this_rand_order = randperm(length(add_cue));
    rand_order = [rand_order; this_rand_order];
    cue_colours = [cue_colours; add_cue(this_rand_order)];
    cue_colours_string = [cue_colours_string; add_strings(this_rand_order)];
    disp(['Colours for run ' num2str(r) ': ' add_strings{this_rand_order(1)} ', ' add_strings{this_rand_order(2)} ', ' add_strings{this_rand_order(3)}])
end

% counter-balance cue shape
if mod(sj_id-1, 4) < 2
    disp('Stimulation is cued by the flipped square!');
    stim_cue = 1;
elseif mod(sj_id-1, 4) > 1
    disp('Stimulation is cued by the upright square!');
    stim_cue = 0;
end

save(['Logs' filesep 'cue-assignment_subject-' num2str(sj_id) '.mat'], 'point_label_assignment', 'cue_colours', 'cue_colours_string', 'stim_cue');


