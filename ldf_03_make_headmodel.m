%% Constructing the Headmodel for EEG forward modelling in FieldTrip

% Segment MRI, make headmodel, make grid
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script takes an (CTF-ALIGNED!) MRI - possibly preprocessed with ldf_01_preprocess_mri.m
% - and computes a 3-element BEM model for EEG or MEG data source reconstruction.
%%
addpath('/home/langdysfun/irichu/Documents'); % Spreadfigures.m for 'segmentation_correction' external script

%%
% NB! check in ldf_01 if bis field correction is ON/OFF
ldf_00_setup
ldf_01_preprocess_mri
ldf_02_bring_mri_to_headspace

%% Before we start we clear all and rerun our setup script

clear all  %# ok
ldf_00_setup
%%
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures

%% Problematic segmentation pps: 02, 06, 07, 08, 09, 11, 12, 17
% incorporated into ldf_01 and 02 now for these pp (no mri reslicing)

% ONLY RUN FOR s11 and s17
%mri_in = ft_read_mri(projpath.mri);
%mri_in.coordsys = 'spm';
%mri_aligned = mri_in;



%% Segment the MRI, get mesh, save mesh (run right up to headmodel code)

if str2num(id) == 17 | str2num(id) == 6 | str2num(id) == 7 | str2num(id) == 8 | ...
        str2num(id) == 9 | str2num(id) == 11  

    segmentation_correction % run another script

else

    % First we load the MRI that was aligned to CTF space:
    load(projpath.mri_aligned);
    % Segment the MRI 
    % This takes time, so be prepared to wait some time. Go read a paper or
    % have a cup of tea :)
    % Remember to RECOMPUTE in case you change the input MRI, e.g. by
    % reslicing or bias field correcting
    
    % If this freaks out with an "Unrecognized field name 'old'." error,
    % then the problem is several spm versions on the path. It might help
    % to restore the defaultpath, otherwise restart MATLAB.
    % Usually this happens if a full older SPM version was added to the path,
    % e.g. for coregistration of data and MRI data via NutMEG.
    
    cfg = [];
    cfg.spmversion = 'spm12';
    cfg.spmmethod = 'old';
    cfg.output = {'brain', 'skull', 'scalp'}; %
    seg = ft_volumesegment(cfg, mri_aligned);
    
    save(projpath.seg, 'seg');
    
    % IGNORE:
    % cfg = [];
    % cfg.spmversion = 'spm12';
    % cfg.spmmethod = 'old';
    % cfg.output = {'gray', 'white', 'csf', 'skull', 'scalp'}; %
    % seg2 = ft_volumesegment(cfg, mri_aligned);
    % 
    % % plot them find which are 1,2,3 (no need, called seg.brain, etc) replace skull and scalp in seg with
    % % skull and scalp from seg2
    % 
    % isequal(seg.skull, seg2.skull) % check if these values are diff
    % isequal(seg.scalp, seg2.scalp)
    % 
    % seg.skull = seg2.skull;  
    % seg.scalp = seg2.scalp;
    
    %% Solving segmentation: deislanding: DROP
    % load(projpath.seg);
    % do sometime before ft_prepare_mesh: deislands3d: remove isolated islands for 3D image (for each slice)
    % for these to run, iso2mesh should be loaded (cfg.method = 'iso2mesh';)
    % can run ft_prepare_mesh lines
    %seg.skull = deislands3d(seg.skull);
    %seg.scalp = deislands3d(seg.scalp);
    %seg.brain = deislands3d(seg.brain); % for s007; 
    
    %% Check segmentation by plotting
    
    % add anatomical information to the segmentation
    %seg.transform = mri_aligned.transform; % Sept 6, 2023: commented out (18/26 pp processed with these lines)
    %seg.anatomy   = mri_aligned.anatomy;
    
    % plot all three tissue types:
    tissue =  {'brain', 'skull', 'scalp'};
    %tissue =  {'scalp'};
    for ii = 1:length(tissue)
    
        cfg = [];
        cfg.method = 'ortho';
        cfg.funparameter = tissue{ii};
        cfg.colorbar = 'no';  % no functional data
        ft_sourceplot(cfg, seg);
        temp = [dirfig 'segmentation_for_mesh_sub-' id '_' tissue{ii} '.png']; print('-dpng',temp);
    end
    
    %% Prepare the mesh
    % We prepare a densely sampled mesh, which we downsample again below.
    cfg = [];
    cfg.tissue = {'brain', 'skull', 'scalp'};
    cfg.method = 'iso2mesh';  
    cfg.spmversion = 'spm12';
    cfg.numvertices = 1e4;   
    mesh = ft_prepare_mesh(cfg, seg);
    
    %% Downsample and repair the meshes
    
    % Check and repair individual meshes using iso2mesh. 
    % This functions are part of the iso2mesh toolbox that ships with FieldTrip. 
    % If MATLAB does not find it, try adding 
    % insert_your_path_to_fieldtrip/external/iso2mesh 
    % to your MATLAB path. 
    % And make sure you are using a reasonably new Fieldtrip version (i.e., 
    % not years old).
    % This is following a routine described here: 
    % https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_mriproc.m
    
    if str2num(id) == 11 % | str2num(id) == 17
        target_sizes = [3000, 3000, 3000];
    else
        target_sizes = [1000, 1000, 1000];
    end
    
    for ii = 1:length(mesh)
        [mesh(ii).pos, mesh(ii).tri] = meshresample(mesh(ii).pos, ...
                                                     mesh(ii).tri, ...
                                                     target_sizes(ii)/...
                                                     size(mesh(1).pos, 1));
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
                                                     mesh(ii).tri, 'dup');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
                                                     mesh(ii).tri, 'isolated');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
                                                     mesh(ii).tri, 'deep');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
                                                     mesh(ii).tri, 'meshfix');
    end
    
    %% Plot the meshes
    
    % Always remember to plot your meshes - make sure there are no holes in any
    % of the meshes by turning the 3D objects. 
    % The different layers get plotted one at a time, press any key to add the
    % next layer to the plot.
    
    figure;
    ft_plot_mesh(mesh(1), 'facecolor', 'red');
    hold on; pause
    ft_plot_mesh(mesh(2), 'facealpha', 0.25, 'facecolor', 'blue');
    hold on; pause          
    ft_plot_mesh(mesh(3), 'facealpha', 0.25, 'edgecolor', [0.8, 0.8, 0.8]);
    view([0 90 0])
    
    temp = [dirfig 'mesh_sub-' id '.png']; print('-dpng',temp);
    
    %% Save mesh
    
    % If this all looks good, we save the mesh - we need it again for the
    % coregistration step
    
    save(projpath.mesh, 'mesh', '-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later

end
%% 
load(projpath.mesh);

%% Make the volume conductor model

% This takes time! -- as in: hours!
tic
cfg = [];
cfg.method = 'dipoli';  % dipoli, bemcp, or openmeeg (the latter needs installation)
cfg.conductivity = [0.33 0.0041 0.33];
vol = ft_prepare_headmodel(cfg, mesh);
toc

save(projpath.vol, 'vol', '-v7.3');

%%
load(projpath.vol);
%% We can double check again if that aligns well with our MRI

if str2num(id) == 2 | str2num(id) == 17 | str2num(id) == 6 | str2num(id) == 7 | str2num(id) == 8 | ...
        str2num(id) == 9 | str2num(id) == 11 | str2num(id) == 12 
    orig_mri = ft_read_mri(projpath.mri); % run ldf_00
    orig_mri.coordsys = 'spm';

    cfg = [];
    cfg.intersectmesh = vol.bnd;
    ft_sourceplot(cfg, orig_mri);

else    
    cfg = [];
    cfg.intersectmesh = vol.bnd;
    ft_sourceplot(cfg, mri_aligned);
end

temp = [dirfig 'vol_aligned_sub-' id '.png']; print('-dpng',temp);


%% Create warped MNI grid

% Let's do the grid warp again ....

% This is what you do if you run a group study. This code takes an
% MNI grid that is shipped with FieldTrip, and warps the positions of this
% template grid to the MRI we supply. That way, we do the source
% reconstruction in individual space (since we warp the grid to the indiv.
% MRI) but they all represent the same coordinates in the MNI brain - which
% makes the positions all comparable to each other in a group analysis.
% One thing you have to keep in mind with this approach is that we cannot
% plot the data on the individual MRI anymore. The data is represented in
% that coordinate space - but the source points do not span a regular grid
% here (since it is warped ...) - and our plotting functions only support
% regular spacing. But: we can just plot on a template brain instead! More
% on that in the plotting file after source reconstruction!

template_grid = load(fullfile(toolboxes.fieldtrip, ...
    '/template/sourcemodel/standard_sourcemodel3d10mm.mat'));

cfg = [];
    if str2num(id) == 2 | str2num(id) == 17 | str2num(id) == 6 | str2num(id) == 7 | str2num(id) == 8 | ...
            str2num(id) == 9 | str2num(id) == 11 | str2num(id) == 12 
    cfg.mri = orig_mri;   % non-aligned that was used
    else
    cfg.mri = mri_aligned;
    end
cfg.grid.template = template_grid.sourcemodel;
cfg.warpmni = 'yes';  % this is unlocking the warping step
cfg.nonlinear = 'yes';  % we warp non-linearly
grid = ft_prepare_sourcemodel(cfg);

save(projpath.grid, 'grid');

