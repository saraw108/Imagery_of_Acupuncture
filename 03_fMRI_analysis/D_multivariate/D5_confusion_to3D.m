% function to save results from the confusion matrices as 3D nifti
% and then probably normalize and smooth them to do tests on them i guess

function D5_confusion_to3D(path_to_all, path_to_confusion, path_to_template, cell_titles)

% in:
% 0. path to all
% 1. confusion matrix
% 2. orgiginal accuracy(_minus_chance) OR mask nifti
% 3. naming of cells as a cellstr, e.g. {'cond1cor', 'cond1ascond2', ...}
% !! naming order: first cycle through rows within one column, then through columns !!
% what it does:
% converts entries (or operations thereof) of confusion matrix to 3D files
% saves all the 3D files

cd(path_to_all);

load(path_to_confusion); % they are name "result"
confusion_cell = results.confusion_matrix.output;
indices = results.mask_index;

template = niftiread(path_to_template);
template_info = niftiinfo(path_to_template);

cond_num = results.n_cond;
temp = zeros([size(template) cond_num^2]);


for l = 1:length(confusion_cell)

    this_voxel = confusion_cell{l,1};
    [x, y, z] = ind2sub(size(template), indices(l,1));

    counter = 0;

    for c = 1:cond_num % these are the COLLUMNS!!

        for n = 1:cond_num

            counter = counter + 1;
            temp(x,y,z,counter) = this_voxel(n,c);

        end

    end

end

for s = 1:cond_num^2

    this_3D = cast(temp(:,:,:,s), 'single');
	disp(['writing output from confusion matrix ' cell_titles{s} '...'])
    niftiwrite(this_3D, [path_to_all filesep 'percent_' cell_titles{s}], template_info);

end

end