%%
cd /project/2422120.01/scripts/source_analysis/Source_recon_langdysfun;
%%
addpath (genpath ('/project/2422120.01/scripts/f_tions/')); %to get the extra f-tions
addpath (genpath ('/home/langdysfun/irichu/Documents/cmocean-main/'));

%%
ldf_00_setup

%%
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures
cmap = colormap(cmocean('thermal'));

%% Load files
load(projpath.fwd) % the forward model
load(projpath.vol) % the headmodel again (to plug into ft_sourceanalysis)
load(projpath.elec_aligned) % load aligned electrodes: use elec_realigned

% get clean EEG data:
sbjnumb = id(2:3); %two last digits from id
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


if any(strcmp(data_reref.label, 'FP1')) %EEG from DCC might erroneously have a channel named "FP1"
    % find where the channel sits
    idx_fp1 = find(strcmp(data_reref.label, 'FP1'));

    % rename it
    data_reref.label{idx_fp1} = 'Fp1';
end

%% CSD matrix
%The beamformer is based on an adaptive spatial filter. For the DICS method, the spatial filter is derived from
% the frequency counterpart of the covariance matrix: the cross-spectral density matrix. 

% CSD is computed from the Fourier transformed data of the single trials. 

% Common spatial filter
%The statistical null hypothesis for both options within (1) is that the data are the same in both conditions,
% and thus the best spatial filter would be the one that is computed using both data conditions together 
% (also known as common filters). This common filter is then applied separately to each condition. 
% To calculate the common filter, we will use the extracted time window, pooled over both conditions.

% Select time window of interest
cfg = [];
cfg.toilim = [-0.8 0];
data_prepic = ft_redefinetrial(cfg, data_reref);

% select frequency range here:
cfg = [];
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.output       = 'powandcsd';
cfg.keeptrials   = 'yes'; % to separate later into conditions
% JM: ft_sourceanalysis does not work with multiple frequency bins, and will probably 
% (silently, and perhaps unexpectedly) average across frequencies in some way.
%cfg.foilim          = [8 30]; 
%cfg.foi = 19; %midpoint between 8 and 30Hz 
cfg.foi = 10; %better choose a point instead of a range
cfg.tapsmofrq    = 4; % amount of smoothing when multitapeting: 4 Hz means plus-minus 4 Hz, i.e. a 8 Hz smoothing box

% for common filter over conditions: use all the data
powcsd_all      = ft_freqanalysis(cfg, data_prepic);


% for conditions:
if task == 1
    cfg.trials = (data_prepic.trialinfo(:,1)==30); % constr verbs 30
elseif task == 2
    cfg.trials = (data_prepic.trialinfo(:,1)==10); % constr noun 10
end
powcsd_con     = ft_freqanalysis(cfg, data_prepic);

if task == 1
    cfg.trials = (data_prepic.trialinfo(:,1)==40); % unconstr verbs 40
elseif task == 2
    cfg.trials = (data_prepic.trialinfo(:,1)==20); % unconstr noun 20
end
powcsd_unc     = ft_freqanalysis(cfg, data_prepic);

%% Compute the spatial common filter and apply it to the conditions
% specify that we would like to keep the computed spatial filter in the output by setting 
% cfg.dics.keepfilter to yes, so that way we can reuse it later.

cfg              = [];
cfg.method       = 'dics';
cfg.sourcemodel  = leadfield;
cfg.headmodel    = vol;
cfg.elec = elec_realigned; 
cfg.frequency    = 10; % single number
cfg.dics.fixedori = 'yes';
cfg.dics.weightnorm = 'unitnoisegain'; % normalise here (e.g. when no control/baseline); not need when you have conditions 
%cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';  % We want to reuse the calculated filter later on
source_all = ft_sourceanalysis(cfg, powcsd_all);

% try running without unitnoisegain
%% Apply common filter to both conditions
%run ft_sourceanalysis again - for both conditions - but specify that we want to use the filter we just computed.
cfg              = [];
cfg.method       = 'dics';
cfg.elec = elec_realigned; 
cfg.sourcemodel  = leadfield;
cfg.headmodel    = vol;
%cfg.frequency    = 10; % don't use this?
cfg.sourcemodel.filter = source_all.avg.filter;  % We apply the previously computed spatial filter

source_con  = ft_sourceanalysis(cfg, powcsd_con);
source_unc = ft_sourceanalysis(cfg, powcsd_unc);

%basically, multiplies data with the filter

%% Compute the difference between conditions
source_diff = source_con;
source_diff.avg.pow = (source_con.avg.pow - source_unc.avg.pow) ./ ...
((source_con.avg.pow + source_unc.avg.pow)/2);
% TFR diff calculated as: '(x1 - x2)/ (x1+x2)'

%% Interpolate onto structural MRI:

% load the template:
template_grid = load(['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
    'template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% load template MRI:
path_template_mri = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        '/template/anatomy/single_subj_T1.nii'];
% read template MRI (MNI space)
temp_mri = ft_read_mri(path_template_mri);
temp_mri.coordsys = 'mni';

% load atlas so that the source is labelled:
    path_aal_atlas = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        'template/atlas/aal/ROI_MNI_V4.nii'];
    atlas = ft_read_atlas(path_aal_atlas);


source_diff.pos = template_grid.sourcemodel.pos;

cfg = [];
cfg.parameter = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int = ft_sourceinterpolate(cfg, source_diff, temp_mri);

% Plot the result
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.funcolorlim   = 'maxabs';
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, source_diff_int);

% Britta's code:
cmo_therm = cmocean('thermal');

    cfg = [];
    cfg.atlas = atlas;
    cfg.funcolormap = cmo_bal; %cmo_therm;
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.funcolorlim = 'maxabs'; %[0, max_pow];  % saturate occipitally
    %cfg.location = source_diff.pos(box_idx(max_idx), :) .* 10; % convert to mm
    cfg.location =  'min';
    ft_sourceplot(cfg, source_diff_int);
    %title(sprintf('sub-%s: %.2f-%.2f sec', sub, cov_time(1), cov_time(2) ))
    %set(gcf, 'color', 'w')


cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
%cfg.funcolorlim    = [0.0 maxval];
cfg.funcolormap    = 'parula'; %cmap; %
%cfg.opacitylim     = [0.0 maxval];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg, source_diff_int);
view ([ 270 0]) % left hemi

%% Grand averaged sources
addpath (genpath ('/project/2422120.01/scripts/f_tions/')); %to get the extra f-tions
%%

sbjlist = {
    '001'
    '002'
    '003'
    '004'
    '005'
    '006'
    '007'
    '008'
    %'09' % leave out for now
    '011'
    '012'
    '013'
    '014' 
    '015'
    '016'
    '017'
    '018'
    '020'
    '021'
    '022'
    '023'
    '024'
    '025'
    '026'  
    '027'
    '028'
 };

task = 1;

% Compute and save sources for GA:
for i = 1:length(sbjlist)

    sbjn = sbjlist(i);
    sbjnumb = cell2mat(sbjn); % convert from cell to mat;does that inside function now
    
    [source_diff] = getSOURCE(task, sbjnumb)

        sources_all{i} = source_diff;

   clear source_diff     
end

if task == 2 
        save(('/project/2422120.01/BIDS/source_noun_all_10Hz'),'sources_all');
elseif task == 1
        save(('/project/2422120.01/BIDS/source_verb_all_10Hz'),'sources_all');
end

%%
load('/project/2422120.01/BIDS/source_noun_all_10Hz.mat')
load('/project/2422120.01/BIDS/source_verb_all_10Hz.mat')
%%

% not to get error: The input sources vary in the field pos
% load the template:
template_grid = load(['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
    'template/sourcemodel/standard_sourcemodel3d10mm.mat']);

for i = 1:length(sbjlist)
    sources_all{i}.pos = template_grid.sourcemodel.pos;
end

% Calculate GAs:
cfg              = [];
%cfg.keepindividual = 'yes'; % default 'no';  store individual results instead (but still gives ga on plots)
cfg.parameter       = 'pow';
ga_source = ft_sourcegrandaverage(cfg, sources_all{:}); % all cells

% calculate GAs on source_difference, then only interpolate here for
% plotting
%% Interpolate onto structural MRI:

% load template MRI:
path_template_mri = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        '/template/anatomy/single_subj_T1.nii'];
% read template MRI (MNI space)
temp_mri = ft_read_mri(path_template_mri);
temp_mri.coordsys = 'mni';


cfg = [];
cfg.parameter = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int = ft_sourceinterpolate(cfg, ga_source, temp_mri);


% Plot GAs:
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.funcolorlim   = 'maxabs';
cfg.opacitymap    = 'rampup';
cfg.funcolormap = 'cool'; % 'hot'
cfg.location = 'min';
ft_sourceplot(cfg, source_diff_int);

% Britta's code for plotting:
cmo_bal = cmocean('balance');
cmo_therm = colormap(cmocean('thermal'));


%[max_pow, max_idx] = max(ga_source.pow);
%[min_pow, min_idx] = min(ga_source.pow);

% load atlas so that the source is labelled:
    path_aal_atlas = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        'template/atlas/aal/ROI_MNI_V4.nii'];
    atlas = ft_read_atlas(path_aal_atlas);

    cfg = [];
    cfg.atlas = atlas;
    cfg.funcolormap = cmo_bal;
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.funcolorlim = 'maxabs'; %[min_pow, 0]; %[0, max_pow];  % saturate occipitally
    %cfg.location = source_diff.pos(box_idx(max_idx), :) .* 10; % convert to mm
    cfg.location =  'min';
    ft_sourceplot(cfg, source_diff_int);



cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = 'maxabs'; %[min_pow max_pow]; %[0.0 maxval];
cfg.funcolormap    = 'cool'; %'cool'; %cmo_therm; % 'parula'
cfg.opacitylim     = 'maxabs'; %[min_pow max_pow]; %[0.0 maxval];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg, source_diff_int);
view ([ 270 0]) % left hemi


%%
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.funcolorlim   = 'maxabs';
cfg.funcolormap    = cmo_bal;
cfg.maskparameter = cfg.funparameter;
%cfg.funcolorlim   = [0.0 1.2];
%cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, ga_source);


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [-0.2 0.2];
cfg.funcolormap    = 'jet'; %cmo_bal; %'jet';
cfg.opacitylim     = [-0.2 0.2];
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfdownsample = 10;  % downsample to speed up processing
ft_sourceplot(cfg, ga_source);
view ([270 0])  

%% Parallel processing
addpath /home/langdysfun/irichu/Documents/toolboxes/fieldtrip/qsub/
%%
% qsubcellfun is a fieldtrip f-tion;
% QSUBCELLFUN applies a function to each element of a cell-array. The
% function execution is done in parallel using the Torque, SGE, PBS or
% SLURM batch queue system.
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = qsubcellfun(fname, x1, x2, 'memreq', 1024^3, 'timreq', 300);


%beware NOT to do something like clear all, clear mex, or clear functions inside your function, as this messes up the environment in which the job gets executed

sbjlist = {
    '001'
    '002'
    '003'
    '004'
    '005'
    '006'
    '007'
    '008'
    %'09' % leave out for now
    '011'
    '012'
    '013'
    '014' 
    '015'
    '016'
    '017'
    '018'
    '020'
    '021'
    '022'
    '023'
    '024'
    '025'
    '026'  
    '027'
    '028'
 };

% task = {'2', '2'}; %don't use

% run here, spawns one extra process per (in this case) subject in subjects:
mem_req = 10 * 1024^3;  % memory requirement in (GB * 1024^3) bytes 
wall_time = 5 * 60;  % time requirement in (minutes * 60) seconds

% @name_of_function 
sources_all_noun = qsubcellfun(@getSOURCEnoun, sbjlist, 'memreq', mem_req, ...
    'timreq', wall_time, 'display', 'yes')

sources_all_verb = qsubcellfun(@getSOURCEverb, sbjlist, 'memreq', mem_req, ...
    'timreq', wall_time, 'display', 'yes')

%%
 save(('/project/2422120.01/BIDS/source_noun_all_10Hz_qsub'),'sources_all_noun');
 save(('/project/2422120.01/BIDS/source_verb_all_10Hz_qsub'),'sources_all_verb');

%%
load('/project/2422120.01/BIDS/source_noun_all_10Hz.mat') % using the loop, size 45.9 mb
load('/project/2422120.01/BIDS/source_verb_all_10Hz.mat')

load('/project/2422120.01/BIDS/source_noun_all_10Hz_qsub.mat') % size 6.4 mb
load('/project/2422120.01/BIDS/source_verb_all_10Hz_qsub.mat')
%%
template_grid = load(['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
    'template/sourcemodel/standard_sourcemodel3d10mm.mat']);

for i = 1:length(sbjlist) %subjects
    sources_all_noun{i, 1}.pos = template_grid.sourcemodel.pos;
end

% Calculate GAs:
cfg              = [];
%cfg.keepindividual = 'yes'; % default 'no';  store individual results instead (but still gives ga on plots)
cfg.parameter       = 'pow';
ga_source = ft_sourcegrandaverage(cfg, sources_all_noun{:}); % all cells

% load template MRI:
path_template_mri = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        '/template/anatomy/single_subj_T1.nii'];
% read template MRI (MNI space)
temp_mri = ft_read_mri(path_template_mri);
temp_mri.coordsys = 'mni';

cfg = [];
cfg.parameter = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int = ft_sourceinterpolate(cfg, ga_source, temp_mri);

% load atlas so that the source is labelled:
    path_aal_atlas = ['/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/' ...
        'template/atlas/aal/ROI_MNI_V4.nii'];
    atlas = ft_read_atlas(path_aal_atlas);

cmo_bal = cmocean('balance');
cmo_therm = colormap(cmocean('thermal'));

cfg = [];
cfg.atlas = atlas;
cfg.funcolormap = cmo_bal;
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs'; %[min_pow, 0]; %[0, max_pow];  % saturate occipitally
%cfg.location = source_diff.pos(box_idx(max_idx), :) .* 10; % convert to mm
cfg.location =  'min';
ft_sourceplot(cfg, source_diff_int);

%% for verb:
for i = 1:length(sbjlist) %subjects
    sources_all_verb{i, 1}.pos = template_grid.sourcemodel.pos;
end

% Calculate GAs:
cfg              = [];
%cfg.keepindividual = 'yes'; % default 'no';  store individual results instead (but still gives ga on plots)
cfg.parameter       = 'pow';
ga_source_v = ft_sourcegrandaverage(cfg, sources_all_verb{:}); % all cells

cfg = [];
cfg.parameter = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int_v = ft_sourceinterpolate(cfg, ga_source_v, temp_mri);

cfg = [];
cfg.atlas = atlas;
cfg.funcolormap = cmo_bal;
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'maxabs'; %[min_pow, 0]; %[0, max_pow];  % saturate occipitally
cfg.location =  'min';
ft_sourceplot(cfg, source_diff_int_v);
