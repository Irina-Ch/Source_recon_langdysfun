%% Forward modelling in FieldTrip

% Segment MRI, make headmodel, make leadfield
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script checks the coregistration and computes the forward model.
%% to get functional data:
addpath (genpath ('/project/2422120.01/scripts/f_tions/')); %to get the extra f-tions

%%
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures

%% Load files

% For the leadfield computation, we need to load the following:
% the electrode information, elec
% the headmodel, vol
% the source grid, grid
% This information will all be combined to compute the forward model, using
% the funciton ft_prepare_leadfield

% read in the headmodel
load(projpath.vol)

% read in the aligned electrodes file
load(projpath.elec_aligned);
% Aug 3 2023: load non-aligned electrodes (to save structure to data_pic_clean before reref)
%load(projpath.elec) % USE ELEC_ALIGNED (the same info there)

% load the source grid
load(projpath.grid)

%%  Sanity check: plot all together (NEVER SKIP THIS STEP!!!)

% Before we continue, we want to make sure that all our forward model
% "ingredients" are in the same coordinate space and well aligned.
% Thus, we plot them together:

figure; hold all;
ft_plot_headmodel(vol, 'edgecolor', 'none', 'facecolor', 'cortex')
alpha 0.5
ft_plot_sens(elec_realigned, 'coilshape', 'point', 'style', 'r.'); % used to be: elec_aligned !
ft_plot_mesh(grid.pos(grid.inside,:));
view([0 90 0])
set(gcf, 'color', [1 1 1])

temp = [dirfig 'forward_model_sub-' id '.png']; print('-dpng',temp);



%%
% all(data_pic_clean.elec.tra == data_reref.elec.tra) % all(): check if all elements are true
% 
% all(elec_aligned.tra == data_reref.elec.tra) % all(): check if all elements are true
% sum(any(data_reref.elec.tra - elec_aligned.tra))
% 
% 
% data_reref.elec2 = rmfield(data_reref.elec, "tra"); % remove tra, see if leadfield still computes
% 


%% Prepare leadfield
% Finally, we can compute the forward model:
% we supply the source grid, the headmodel (vol) and the electrode
% information (elec).
% We tell the algorithm to not reduce rank (this is EEG data, so we are
% able to see radial and tangential sources, thus we use the full rank of
% the leadfield) and to not normalize the leadfield (if we do want to get
% rid of the center-of-head bias, we prefer to normalize the weights of the
% beamformer:Westner et al. 2022. DOI: 10.1016/j.neuroimage.2021.118789 for more info).

cfg = [];
cfg.sourcemodel = grid;
cfg.headmodel = vol;
cfg.elec = elec_realigned; %elec_aligned; 
cfg.channel = 'EEG'; %
cfg.reducerank = 'no';
cfg.normalize = 'no';
leadfield = ft_prepare_leadfield(cfg);

save(projpath.fwd, 'leadfield')
