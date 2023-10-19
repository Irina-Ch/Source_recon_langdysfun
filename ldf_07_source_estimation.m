%% Source estimation

% Segment MRI, make headmodel, make leadfield
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script shows an example source reconstruction using a beamformer.
% So far, the scripts have been very universal and can be used for many
% different source reconstruction scenarios. This script is where the data
% comes into play, so you will have to make your own decisions on different
% parameters. 
% For demonstration purposes, I assume the following:
% We recorded a dataset with visual stimulation and we want to just see the
% activity after the stimulus was presented.
%%
addpath (genpath ('/project/2422120.01/scripts/f_tions/')); %to get the extra f-tions
addpath (genpath ('/home/langdysfun/irichu/Documents/cmocean-main/'));

%%
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures

%% Load files

% the forward model
load(projpath.fwd)

% we also need to supply the headmodel again:
load(projpath.vol)

% load aligned electrodes
load(projpath.elec_aligned)

%% Load functional data
% We want to run the forward model on rereferenced-to-all data too, so:
% 1. load functional, save original 'elec' structure there, reref it, use
% data_reref.elec from now on!

%NB!! REMEMBER TO CHANGE SUBJ NUMBER
% It is assumed the data is cleaned and epoched around stimulus onset - or
% that you have a script that can produce such data. ORIGINAL SCRIPT: struct 'epochs'

%str2num(sbjnumb) == str2num(id) 
sbjnumb = '28';
task =  2; % 1 for verb, 2 for noun
layout = '/project/2422120.01/BIDS/ActiCap_64Ch_DCC_customized.mat';
[data_clean, badtrials] = get_data_clean(task, sbjnumb); 

cfg=[];
cfg.artfctdef.visual.artifact = badtrials;
data_pic_clean = ft_rejectartifact(cfg, data_clean); % data_pic if segmented to sentence onset

cfg=[];
cfg.implicitref   = 'TP9'; % online ref, left 
cfg.reref         = 'yes';
cfg.refchannel    = {'all'}; %{'TP9' 'TP10'}; % cell-array with new EEG reference channel(s), 'all' for a common average reference
cfg.refmethod     = 'avg'; % default; other options: 'median', 'rest', 'bipolar' or 'laplace' 
data_reref        = ft_preprocessing(cfg, data_pic_clean); 


% NOTE:
% Some files from the DCC might erroneously have a channel named "FP1",
% while it actually should be named "Fp1". If this is the case, we need to
% fix this here (or even earlier in your pipeline as this will also create
% mismatches for topo-plotting).

if any(strcmp(data_reref.label, 'FP1'))
    % find where the channel sits
    idx_fp1 = find(strcmp(data_reref.label, 'FP1'));

    % rename it
    data_reref.label{idx_fp1} = 'Fp1';
end


% plot VEP too:
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;
cfg.lpfilttype = 'firws';
datafilt = ft_preprocessing(cfg, data_reref);

cfg=[];
cfg.channel = [1:64]; % exclude external electrodes
tlk.vep = ft_timelockanalysis(cfg, datafilt);

cfg = [];
cfg.layout = layout;
cfg.colormap = cmocean('balance');  
cfg.channel = {'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' ...
    'P5' 'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'PO7' 'PO9' 'PO8' ...
    'PO10' 'PO3' 'PO4' 'POz' 'O1' 'Oz' 'O2'};
cfg.xlim   = [-0.15 0.25]; % plot -150 to 250 ms
cfg.ylim = 'maxabs';
ft_singleplotER(cfg, tlk.vep); % click on peak 50-100ms to plot (visual areas)

cfg = []; % plot VEP: time-window of around 100ms
cfg.colormap = cmocean('balance');  
cfg.layout = layout;
cfg.ylim = 'maxabs';
cfg.xlim   = [0.085 0.125];  % [0 : 0.05 : 0.3] as a vector to plot a series of topos
ft_topoplotER(cfg,  tlk.vep); colorbar

% IMPORTANT:
% Remember that this data needs to be average-referenced! YES

% We also remove the EOG and EMG channels from the data: ALREADY DONE
% cfg = [];
% cfg.channel = {'EEG', '-HEOG', '-VEOG', '-EMG'};
% epochs = ft_selectdata(cfg, epochs);

% NB!!! NOT ENOUGH JUST TO RENAME/REMOVE
% REMOVED TP9 which was returned by rereferencing from Ref (making it 65 in data_reref channels instead of 64 in the leadfields )
%setdiff(data_reref.label, leadfield.label) % values in data that are not in leadfield: 'TP9'
%strfind(leadfield.label,'TP9')
% cfg = [];
% cfg.channel = {'EEG', '-TP9'};
% data_reref = ft_selectdata(cfg, data_reref);

% Doing this instead:
% 1. data_pic_clean.elec = elec; (before reref, save the original elec structure)
% 2. re-reference that
% 3. cfg.elec = data_reref.elec (for ft_sourceanalysis)

%% Compute the covariance matrix for computing the beamformer

% For your own data, you will likely have to adjust the covariance window
cfg = [];
cfg.covariancewindow = [-0.15 0.15]; % ORIGINAL: common filter for visual time window!  
%Common filter has activity from both conditions (or both windows in this case)
%cfg.covariancewindow = [0 0.2]; % to compute visual activity (instead of common filter; try diff time windows depending on pps VEP)
cfg.covariance = 'yes';
cov = ft_timelockanalysis(cfg, data_reref);

% we compute the rank of the covariance matrix:
cov_rank = rank(cov.cov)

%% Septermber 8, 2023: check cov matrix and leadfields CODE DOESNT WORK
figure; imagesc(cov.cov); colorbar; 
caxis([-max(max(cov.cov)) max(max(cov.cov))]); % so that it is centered around 0

%coords 18, -98, 6 (right calcarine)
% to plot leadfields
leadfield
%dim contains the dimensions of the grid in which the 6804 (18x21x18) sources are placed
%pos contains the xyz coordinates for the sources in the source model
% inside contains a logical vector indicating whether the source is a source or not
%leadfield a cell array, where each cell (6804) contains a matrix of 65 rows (n channels) 
% and 3 columns (_xyz-coordinates). (The ones that are outside the brain are empty though). 
%leadfielddimord indicates the dimension. Each cell is a position, which contains a matrix ordered
% by channels x orientations (xyz)


size(leadfield.pos) %6804x3 position (orientation)
size(leadfield.leadfield) % 1 x 6804


leadfield.pos(1,1) % position of voxels, find the closest to 18, -98, 6

[val1, indx1] = min(abs(leadfield.pos(:,1)-18)) % find nearest number
leadfield.pos(indx1) %17.9  

%[val, indx] = min(abs(leadfield.pos(:)-ourcoord)) % find nearest number
%leadfield.pos(indx)

ourcoord = [18, -98, 6];
[k, dist] = dsearchn(leadfield.pos, ourcoord) % dsearchn() nearest point search

closest = leadfield.pos(k,:) % closest?  18.3093  -78.2464    1.7057

% leadfield.leadfield{1,2000}(1,:) % cell array {1, cell number}(indices within that cell)

leadfield.leadfield{k} % 1008 this is individual position, not in mni?

ourcoord = [18, -98, 6];
%[k, dist] = dsearchn(template_grid.sourcemodel.pos, ourcoord) % in mni; % use the index you get on the source on the leadfield 
%closest = template_grid.sourcemodel.pos(k,:) % ?? 8.5000  -12.0000    6.0000
%leadfield.leadfield{k} 

%inside_sources = find(leadfield.inside);
pos_mm = ft_convert_units(template_grid.sourcemodel, 'mm') % convert to mm since leadfield is in mm
[k, dist] = dsearchn(pos_mm.pos, ourcoord) 
closest = pos_mm.pos(k,:) % 15  -100    10
size(leadfield.leadfield{k} )
leadfield.leadfield{k}(:,1) % x; col2: y; col3: z

data = [];
data.dimord = 'chan_time';
data.label = data_reref.label;
%data.trialinfo = 1;
data.time = 1;

data.avg = leadfield.leadfield{k}(:,1); % X
cfg = [];
cfg.colormap = cmocean('balance');  
cfg.layout = layout;
%cfg.channel = 'all';
cfg.zlim = 'maxabs'; 
cfg.style = 'straight'; %no lines
figure; ft_topoplotER(cfg, data); colorbar;

data.avg = leadfield.leadfield{k}(:,2); % X
cfg = [];
cfg.colormap = cmocean('balance');  
cfg.layout = layout;
%cfg.channel = 'all';
cfg.zlim = 'maxabs'; 
cfg.style = 'straight';
figure; ft_topoplotER(cfg, data); colorbar;

data.avg = leadfield.leadfield{k}(:,3); % X
cfg = [];
cfg.colormap = cmocean('balance');  
cfg.layout = layout;
%cfg.channel = 'all';
cfg.zlim = 'maxabs'; 
cfg.style = 'straight';
figure; ft_topoplotER(cfg, data); colorbar;

X = leadfield.leadfield{k}(:,1)
Y = leadfield.leadfield{k}(:,2)
Z = leadfield.leadfield{k}(:,3)

pos_inside_ind = find(pos_mm.inside);
pos_inside = pos_mm.pos(pos_inside_ind);

pos_mm.pos(template_grid.sourcemodel.inside, 1)

figure; 
%plot3(pos_mm.pos(:,1), pos_mm.pos(:,2),pos_mm.pos(:,3), 'o')
plot3(pos_mm.pos(template_grid.sourcemodel.inside, 1), pos_mm.pos(template_grid.sourcemodel.inside,2), ...
    pos_mm.pos(template_grid.sourcemodel.inside, 3), 'o')
grid on
hold on
%figure; 
%plot3(X,Y,Z, 'o', 'Color', 'r');
plot3(closest(1),closest(2), closest(3), 'o', 'Color', 'r');
grid on
%p.Color = "red";




%% Now we can make the spatial filter

% We compute an LCMV beamformer - that is a beamformer for time domain
% data. 
% We supply the headmodel and the forward model, as well as the electrode
% information.
% OPTIONS:
% We use a weight-normalized beamformer (unit-noise-gain beamformer) and
% compute the dipole orientation (per source point) that maxmimizes output
% power. We handle rank-deficient covariance matrices by using a truncated
% pseudo-inverse (we set kappa to our rank).
% Lastly, we tell the algorithm to save the spatial filter in the output -
% as we want to apply it to different data.

% Referencing Westner et al. 2022, DOI: 10.1016/j.neuroimage.2021.118789 : 
% You can read more about weight normalization in section 2.2, "Weight
% normalization strategies". 
% You can find information about the dipole orienation in section 2.1, "Basic
% beamformer formulations".
% The inversion of rank deficient covariance matrices is discusse in detail
% in 3.1, "Estimation and inversion of covariance matrices".
% Common spatial filters are discussed in 3.7, "Common choices for beamformer
% analysis pipelines".

setdiff(cov.label, leadfield.label) % check if (error: cannot match the channels in the sourcemodel to those in the data)

% MAKE LCMV (for the whole cov matrix, i.e. the common spatial filter)
cfg = [];
cfg.method = 'lcmv';
cfg.headmodel = vol;  % headmodel
cfg.sourcemodel = leadfield;  % forward model
cfg.elec = elec_realigned;  % electrode information (error: cannot match the channels in the sourcemodel to those in the data)
cfg.lcmv.weightnorm = 'unitnoisegain';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.kappa = cov_rank;
cfg.lcmv.keepfilter = 'yes';
spat_filter = ft_sourceanalysis(cfg, cov);

%% Apply spatial filters

% FieldTrip does not have seperate calls for creating and applying a
% spatial filter. The variable spat_filter above in fact contains the
% source reconstruction of the data in cov.
% However, we want to use the spatial filter we created using data from the
% baseline and active window of our data (cfg. the covariance window) and
% apply it to the baseline and active time points *separately* so that we
% can contrast them.
% This approach is also called using a "common spatial filter", since our
% beamformer is not biased toward either the active or baseline activity.

% For this to work, we will have to recompute covariance matrices for both
% the baseline and the active window and the apply the common spatial
% filter to these covariance matrices.

% Compute covaraince matrices for pre- and post-stim: RUN THIS!
cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.15 -0.05]; % [-0.2 -0.05] can try 0.15 instead of 0.2
cov_prestim = ft_timelockanalysis(cfg, data_reref);

cfg.covariancewindow = [0.05 0.15];
cov_poststim = ft_timelockanalysis(cfg, data_reref);

if 0 % DOES NOT RUN code in the statement
% In our case, diff between conditions:
cfg = []; 
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.8 0];
if task == 1
    cfg.trials = (data_reref.trialinfo(:,1)==30); % constr verbs 30
elseif task == 2
    cfg.trials = (data_reref.trialinfo(:,1)==10); % constr noun 10
end
cov_constr = ft_timelockanalysis(cfg, data_reref); % CONSTR


if task == 1
    cfg.trials = (data_reref.trialinfo(:,1)==40); % unconstr verbs 40
elseif task == 2
    cfg.trials = (data_reref.trialinfo(:,1)==20); % unconstr noun 20
end
cov_unconstr = ft_timelockanalysis(cfg, data_reref); % UNCONSTR

end

% APPLY LCMV (multiplies cov matrix by weights)
% apply pre-computed spatial filter to this:
cfg = [];
cfg.method = 'lcmv';
cfg.elec = elec_realigned; % used to be 'elec' 
cfg.headmodel = vol; % not used by ftion (relic of prev version)
cfg.sourcemodel = leadfield; % not used by ftion
cfg.sourcemodel.filter = spat_filter.avg.filter;
source_prestim = ft_sourceanalysis(cfg, cov_prestim);
source_poststim = ft_sourceanalysis(cfg, cov_poststim);

%% Contrast conditions

% now we obtained the source estimates for pre- and post-stim activity.
% Let's contrast them!

% Tip: you can also do this using ft_math()
source_diff = source_prestim;  % make a copy of the source variable
source_diff.avg.pow = (source_poststim.avg.pow - source_prestim.avg.pow); % constrained - unconstrained

%% Prepare plotting

% We are almost there! Just the plotting step is missing.
% First, we have to remember that we are working with a warped grid. Go
% back to the end of ldf_04_make_headmodel.m if you need a refresher of
% what that means.
% Thus, we cannot plot the data on the individual MRI, but have to use the
% template MRI in MNI space!
% Let's load it:

template_mri = fullfile(toolboxes.fieldtrip, ...
    '/template/anatomy/single_subj_T1.nii');

temp_mri = ft_read_mri(template_mri);
temp_mri.coordsys = 'mni';

% Now remember that the positions of our source reconstruction are still in 
% individual space (because our headmodel and electrde positions etc are all in
% that individual space, the resulting source estimation is too) - but they
% correspond 1:1 to coordinates in the MNI space (that is what we have
% carefully set up during the warping).
% So we can go ahead, load the template and replace the source grid
% coordinates! (Yep, sounds crazy but works fine)

template_grid = load(fullfile(toolboxes.fieldtrip, ...
    '/template/sourcemodel/standard_sourcemodel3d10mm.mat'));
% Side note: re-loading things you used earlier in the pipeline is when I
% would strongly recommend to use some kind of setup-file and a path tree
% (like our projpath object). That way you minimize the danger of loading a
% *different* file than before (and then chaos ensues ...).
% Here I kept the full path in the scripts to make it very clear to you
% which file we are loading and from where, but for a "real" analysis I
% would move this info to ldf_00_setup.m to keep it constant across
% scripts.

% Just replace the positions!
source_diff.pos = template_grid.sourcemodel.pos;

%% Interpolate the source estimate

% The MRI has a different resolution than our source estimate, so we need
% to interpolate:
cfg = [];
cfg.parameter = 'avg.pow';
cfg.interpmethod = 'nearest';
source_int = ft_sourceinterpolate(cfg, source_diff, temp_mri);

% And lastly, plot it!

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs';
cfg.locationcoordinates = 'voxel'; % for crosshairs
cfg.location = [55 15 40];% so crosshair starts there: roughly around right calcarine sulcus
ft_sourceplot(cfg, source_int);

%%
if task == 1
    temp = [dirfig 'source_plot_visual_verb_sub-' id '.png']; print('-dpng',temp);
elseif task == 2
    temp = [dirfig 'source_plot_noun_visual_sub-' id '.png']; print('-dpng',temp);
end



