%% Coregister MRI and electrodes using surface matching

% Align MRI and structural scan using surface matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% CREDIT: based on code by Sarang S. Dalal
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ldf_00_setup;
%% Add some necessary extra paths

% We need NutMEG on our path for several functions:
addpath(toolboxes.nutmeg);
addpath(fullfile(toolboxes.nutmeg, 'external'));

%% Load the structural scan
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures

% we read the structural scan as a "headshape"
projpath.struct_scan = ['/project/2422120.01/3d_scans/sub-', id, '/Model_', id, '.obj']; % NB! models are outside BIDS since we wont publish them
head_surface = ft_read_headshape(projpath.struct_scan);

% make sure things are in millimeters and check by plotting:
head_surface = ft_convert_units(head_surface, 'mm'); 
figure;
ft_plot_mesh(head_surface); 

% copy the mesh over:
scan_mesh.vertices = head_surface.pos;
scan_mesh.faces = head_surface.tri;

%% Load the MRI mesh

load(projpath.mesh);

% copy the scalp mesh over: -- NOTE: this is the right index for the scalp
% mesh if you followed the rest of the pipeline. If you got your mesh from
% somewhere else, you might have to adjust the index.
mri_mesh.vertices = mesh(3).pos;
mri_mesh.faces = mesh(3).tri;

% save some memory:
clear head_surface 

%% First, we roughly align the two models

% For the iterative process to work, we need to roughly align the
% structural scan and the MRI mesh. For that, please identify the points
% the figure prompts you to select. You click on the point and press ENTER.
% This is not the coregistration step yet - that one will be done fully
% automated.

[mri_fids, surface_fids] = identify_points(mri_mesh, scan_mesh);

% define initial rigid transform based on identified fiducials
tfm_0 = rigid_coreg(surface_fids, mri_fids);

% apply this transform on the structural scan
scan_mesh.vertices_0 = nut_coordtfm(scan_mesh.vertices, tfm_0);

%% Next, we select the face for automatic coregistration

% We want to use the face only for automatic coregistration: the top and 
% back of the head will contain the electrodes for the structural scan and 
% thus not perfectly align with the MRI head surface (i.e., this would
% introduce coregistration error).
% If the MRI has signal cancellations (e.g. because of metal in the mouth),
% select a point above that region when prompted to select the chin or
% above the mouth. 

face_points = select_face(scan_mesh);

%% Run the automatic coregistration via ICP

% Coregister the two meshes via the iterative closest point algorithm. Use
% only the face points for that.
[tfm, scan_mesh] = make_icp_coreg(mri_mesh, scan_mesh, face_points, tfm_0);

temp = [dirfig 'coreg_face_sub-' id '.png']; print('-dpng',temp);

%% Visualize the coregistration agreement

% We plot the coregistration result. The MRI in blue and the structural
% scan in red. Both are plotting semi-transparent, so you should be able to
% make out both.
% The will not perfectly 100% align (that is not possible), but the should
% clearly share the same space and orientation, sometimes the MRI
% overlapping the structural scan and sometimes the other way around.

figure;

% view angle 1
subplot(1, 2, 1)
trisurf(mri_mesh.faces, ...
    mri_mesh.vertices(:, 1), ...
    mri_mesh.vertices(:, 2), ...
    mri_mesh.vertices(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'blue', 'FaceAlpha', 0.65)
hold on
trisurf(scan_mesh.faces, ...
    scan_mesh.vertices_1(:, 1), ...
    scan_mesh.vertices_1(:, 2), ...
    scan_mesh.vertices_1(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.65)
view([90 0])
axis equal
lightangle(240,45);
lightangle(-240,45);

% view angle 2
subplot(1, 2, 2)
trisurf(mri_mesh.faces, ...
    mri_mesh.vertices(:, 1), ...
    mri_mesh.vertices(:, 2), ...
    mri_mesh.vertices(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'blue', 'FaceAlpha', 0.65)
hold on
trisurf(scan_mesh.faces, ...
    scan_mesh.vertices_1(:, 1), ...
    scan_mesh.vertices_1(:, 2), ...
    scan_mesh.vertices_1(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.65)
view([10 0])
axis equal
lightangle(240,45);
lightangle(-240,45);

temp = [dirfig 'coreg_agreement_sub-' id '.png']; print('-dpng',temp);


%% Apply the coregistration to the electrodes

load(projpath.elec);

if ~strcmp(elec.unit, 'mm')
    elec = ft_convert_units(elec, 'mm');
end

elec_aligned = elec;
elec_aligned.elecpos = nut_coordtfm(elec.elecpos, tfm);
elec_aligned.chanpos = nut_coordtfm(elec.chanpos, tfm);

% re-label Ref to TP9 in elec_aligned
refindx = find(strcmp(elec_aligned.label, 'Ref'));
elec_aligned.label(refindx) = {'TP9'};

%% Visualize the electrodes together with the MRI

figure;
ft_plot_mesh(mesh(3));
hold on
ft_plot_sens(elec_aligned)

temp = [dirfig 'elec_with_MRI_sub-' id '.png']; print('-dpng',temp);

%% Double check the labelling as well

figure;
ft_plot_headshape(mesh(3))
ft_plot_sens(elec_aligned, 'label', 'on', 'fontsize', 15, 'elecshape', 'disc', 'elecsize', 10)

%% Save the aligned electrodes

%elec_aligned.label{1}= 'TP9'; % change Ref to TP9 since it is returned in
% functional data after reref to all, ending in 65 elecs in total (need 65
% here too so that leadfield.label == data_reref.label
%NB! elec.tra also has info about referencing scheme!! so do instead:
% data.elec = elec
% then: reference that
% then: cfg.elec = data_reref.elec

save(projpath.elec_aligned, 'elec_aligned');



%% if elecs are not aligned well + change Ref to TP9
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures

load(projpath.elec_aligned);
load(projpath.mesh);

% re-label Ref to TP9 in elec_aligned
refindx = find(strcmp(elec_aligned.label, 'Ref'));
elec_aligned.label(refindx) = {'TP9'};
%%

cfg           = [];
cfg.method    = 'project'; %'interactive'
cfg.elec      = elec_aligned;
cfg.headshape = mesh(3); % mesh(3) = scalp!
elec_realigned  = ft_electroderealign(cfg);

figure;
ft_plot_mesh(mesh(3));
hold on
ft_plot_sens(elec_aligned);

figure;
ft_plot_mesh(mesh(3));
hold on
ft_plot_sens(elec_realigned);
temp = [dirfig 'elec_with_MRI_REALIGN_sub-' id '.png']; print('-dpng',temp);
% DONT FORGET TO SAVE!
%%

save(projpath.elec_aligned, 'elec_aligned', 'elec_realigned');

%% 
load(projpath.grid);
load(projpath.vol);

figure;
hold all;
ft_plot_headmodel(vol, 'edgecolor', 'none', 'facecolor', 'cortex')
alpha 0.5
ft_plot_sens(elec_realigned, 'coilshape', 'point', 'style', 'r.');
ft_plot_mesh(grid.pos(grid.inside,:));
view([0 90 0])
set(gcf, 'color', [1 1 1])

