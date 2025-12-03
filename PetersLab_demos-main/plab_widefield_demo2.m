%% Widefield demo 2
%
% A demo on how to work with widefield data in batch
% (do plab_widefield_demo first)
%
% This code shows how to align widefield data across animals and mice. Note
% that these functions are likely to change/improve over time.
%
% If you would like to use these functions for regular analysis, I would
% suggest that you create a script that has these functions in a pipeline
% to run for your animals.

%% Github repositories to download

% These repositories are necessary for this demo:
%
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab
% https://github.com/kwikteam/npy-matlab

%% Location for saved alignments

% Widefield alignments are saved as a tranform matrix for each animal and
% day. When an alignment is applied (e.g. when loading data), it loads the
% alignment for that animal and applies the transform matrix for that day. 

% Alignments for all animals are kept in: 
% \\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\widefield_alignment
% 
% This has two folders: 
% animal_alignment - the transform matrices for each animal
% retinotopy - retinotopy for each animal used in across-animal aligning
% 
% There are also general files which are not editable, e.g. master_vfs.mat
% is the 'master' retinotopy that all animals are aligned to

% Keeping all alignments in one place makes the pipeline easier, but it
% also means you need to be careful not to mess up other people's data.
% Hopefully this is unlikely, since everyone will only create/use
% alignments for their own animals.


%% Widefield alignment function

% I use the function 'plab.wf.wf_align' to create and apply widefield
% alignments. Check out the help for 'plab.wf.wf_align' for details about how it
% works. Below are examples for how to use it. 

% The general syntax is: 
% im_aligned = wf_align(im_unaligned,animal,day,align_type,master_align)
% 
% This applies alignments to a set of images 'im_unaligned' to create
% aligned images 'im_aligned', or creates and saves new alignments.


%% Creating alignments 

% We can use this animal as an example (data can be overwritten):
animal = 'AP013';

% There are two steps to creating alignments: 

%------- 1) aligning animal across days

% Images are aligned across days using the average blue image. Aligment
% here is rigid, meaning only translation and rotation, since the scaling
% should not change across days. All images are aligned to the first day,
% making the first day the 'master' image for each animal.

% The function automatically finds and loads these, so no unaligned images
% need to be entered (the first argument is blank). This also grabs the
% date for each image, so no dates need to be entered (the third argument
% is blank). Only the second argument (animal), and third argument
% (alignment type) need to be entered. 

% Create a new across-day alignment for this animal by running the command
% below. This will display "Rigid aligning images..." in the command line
% while this works, and then display a scrollable set of images:
ap.wf_align([],animal,[],'new_days');

% The command window should now display "Manual align? (#,s=save, q=quit)".
% entering 's' will save the alignments as-is, entering 'q' will quit
% without saving, and entering a number (e.g. 3) will let you manually
% align that day. 

% Try scrolling through the images, and note that there is a craniotomy on
% day 2 and 3, so they look different from day 1. Days 1-2 should look
% pretty well aligned, but 3 is slightly off. In the command window, type
% '3' and enter.

% This brings up the control point selection tool. Manually click a few
% points on the left (day 3) and right (day 1 = master) images that match
% each other. For the best alignment, select 4+ points, and make the points
% as widely disributed as possible. When you're finished selecting, close
% the window by pressing 'x'. 

% The aligned images will now be displayed again. If you're happy with this
% alignment, you can enter 's'. If you want to do more manual alignment,
% you can enter the number of the day you want to align.

% Once you enter 's', the command window will display "Saved day transforms
% for AP013."


%------- 2) aligning animal to a 'master' alignment (across animals)

% After across-day alignments have been created, you can now align an
% animal to the 'master' alignment, which lets multiple animals be aligned
% to each other. 

% We do this using visual sign maps (vfs), where we align each animal's vfs
% to a 'master' vfs stored in the alignment folder, which was created from
% the average of many aligned mice. Here is what the master vfs looks like:
load('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\widefield_alignment\master_vfs.mat');
figure;
imagesc(master_vfs);
axis image;
title('Master VFS');

% In widefield demo 1, we used 'plab.wf.retinotopy_vfs' to create the
% visual field sign map. As a reminder, here we can find sparse noise
% recordings for this animal, and create and plot the vfs. There are more
% than one sparse noise recordings, so this will only plot the first:
recordings = plab.find_recordings(animal,[],'sparse_noise');
plot_recording = 1;
vfs = plab.wf.retinotopy_vfs(animal, ...
    recordings(plot_recording).day,recordings(plot_recording).recording{1});

% To automatically find and save the vfs for all sparse noise recordings
% for a given animal, use the function 'plab.wf.retinotopy_vfs_batch',
% which runs the above retinotopy function and saves into the alignment
% path:
plab.wf.retinotopy_vfs_batch(animal);

% When it's finished, the command window should print the location where it
% saved the vfs maps. 

% After the vfs maps are saved, we can create the across-animal alignment.
% This first creates an average vfs by loading in saved vfs maps, aligning
% across days, and averaging. It then aligns the average vfs for that
% animal with the master vfs. Create that alignment with this command:
plab.wf.wf_align([],animal,[],'new_animal');

% When it's finished running, it will display the master and aligned
% images, and the command window will print "Save animal alignment for
% AP013? (y/n):". If you are happy with the alignment, enter 'y'. If you
% are not, you can exit without saving with 'n'. This doesn't have a manual
% option at the moment, so if this alignment isn't good, then it can't be
% corrected. 

% Enter 'y', and the command window should display "Saved new animal
% alignment for AP013."

% After this step, this animal (AP013) will now have both across-day and
% across-animal alignments.

%% Alignment summary

% In summary, your pipeline for alignment should look like this: 

% Define animal
animal = 'AP016';
% Create across-day alignments
plab.wf.wf_align([],animal,[],'new_days');
% Get and save VFS maps for animal
plab.wf.retinotopy_vfs_batch(animal);
% Create across-animal alignments
plab.wf.wf_align([],animal,[],'new_animal');


%% Applying alignments

% Once an alignment is created, it can be applied with the same function by
% entering the recording day. 

% For example, we can load in an average image from AP013 and align it: 
animal = 'AP016';
rec_day = '2024-02-20';
avg_im_filename = plab.locations.filename('server',animal,rec_day,[],'widefield','meanImage_blue.npy');
avg_im = readNPY(avg_im_filename);
avg_im_aligned = plab.wf.wf_align(avg_im,animal,rec_day);

figure;
subplot(1,2,1);
imagesc(avg_im); axis image;
title('Unaligned image');
subplot(1,2,2);
imagesc(avg_im_aligned); axis image;
title('Aligned image');

% In order to align SVD images, we can align the U's (spatial components)
% in the same way, i.e. 'U_aligned =
% plab.wf.wf_align(U_unaligned,animal,rec_day)'. This is automatically done
% in 'ap.load_experiment' if an alignment exists. 

% Once a widefield image is aligned, you can draw cortical areas on top
% using the function 'ap.wf_draw'. This shows the aligned average image,
% and draws the cortical areas on top: 
figure;
imagesc(avg_im_aligned); axis image;
ap.wf_draw('ccf','r');
title('Aligned image with cortical areas');



