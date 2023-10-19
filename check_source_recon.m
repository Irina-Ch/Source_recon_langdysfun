%% Check source recon CTXTVerb
%%
addpath (genpath ('/project/2422120.01/scripts/f_tions/')); %to get the extra f-tions
addpath('/home/langdysfun/irichu/Documents'); % Spreadfigures.m for 'segmentation_correction' external script

%%
% Fieldtrip needs to be on path.
codetree.cmocean = fullfile('/home/langdysfun/irichu/Documents/cmocean-main');
codetree.fieldtrip = fullfile('/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/');
addpath(codetree.cmocean)

%% load general stuff

% load the template
template_grid = load(fullfile(codetree.fieldtrip, ...
    '/template/sourcemodel/standard_sourcemodel3d10mm.mat'));

%% Run across subjects

%subjects = arrayfun(@(X) sprintf('%.3d', X), 1:28, 'uniformoutput', 0)';
%subjects = setdiff(subjects, {'006', '009', '010', '011', '017', '019'}); % remove these

%for s_idx = 1:length(subjects) %Ir: added 1:
%    sub = subjects{s_idx}; % enable the loop to run multiple pps
    
    %use if running 1 pp:
    sub = '006';

    bids_path_anat = ['/project/2422120.01/BIDS/sub-', sub, '/anat'];
    %bids_path_data = ['/project/2422120.01/BIDS/sub-', sub, '/eeg/raw']; %for raw data

    % load headmodelling stuff
    load(fullfile(bids_path_anat, [sub, '_fwd.mat']));  % leadfield
    load(fullfile(bids_path_anat, [sub, '_elec_aligned.mat']));  % elec_aligned
    load(fullfile(bids_path_anat, [sub, '_vol.mat']));  % this loads a lot for 001

%% original code using raw data:
%     % raw data path
%     eeg_path = fullfile(bids_path_data, ['sub-', sub, '_task-noun_eeg.eeg']);  % data
% 
% 
%     % Read in the data:
%     cfg = [];
%     cfg.dataset = eeg_path;
%     cfg.trialdef.eventtype = 'Stimulus';
%     cfg.trialdef.eventvalue   = {'S 10' 'S 20'};
%     cfg.trialdef.prestim = 10;
%     cfg.trialdef.poststim = 1;
%     cfg = ft_definetrial(cfg);
% 
%     cfg.lpfilter = 'yes';
%     cfg.lpfreq = 30;  % for this check
%     cfg.lpfilttype = 'firws';
%     cfg.demean = 'yes';
%     epochs = ft_preprocessing(cfg);
% 
%     % Get rid of non-EEG electrodes
% 
%     cfg = [];
%     cfg.channel = 'EEG';
%     epochs = ft_selectdata(cfg, epochs);
% 
%     % Re-reference
%     cfg = [];
%     cfg.reref = 'yes';
%     cfg.refchannel = 'all';
%     cfg.implicitref = 'TP9';
%     epochs = ft_preprocessing(cfg, epochs);

    sbjnumb = sub(2:3); % take 2 last digits from sub; sbjnumb = '01';
    task =  1; % 1 for verb, 2 for noun
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
    
    
    if any(strcmp(data_reref.label, 'FP1'))
        % find where the channel sits
        idx_fp1 = find(strcmp(data_reref.label, 'FP1'));
    
        % rename it
        data_reref.label{idx_fp1} = 'Fp1';
    end

    %% make cov, check topo and ERF
    cov_time =  [0.05, 0.15]; % no subtraction later

    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = cov_time;
    cov = ft_timelockanalysis(cfg, data_reref); %originally: (cfg, epochs)

    f1 = figure;

    % call once to inialize the y-axis already
    cfg = [];
    cfg.figure = f1;
    cfg.channel = '*O*';  % all occipital channels
    cfg.xlim = [-0.3, 0.5];
    ft_singleplotER(cfg, cov);
    hold on

    % add shaded area for cov time
    ylims = get(gca, 'ylim');
    ar = area(cov_time([1, 1, 2, 2]), [ylims, ylims([2, 1])], ...
        'handlevisibility', 'off');
    ar.FaceColor = [0.95, 0.95, 0.95];
    ar.EdgeColor = 'none';
    ar.BaseValue = ylims(1);
    hold on;

    % re-plot on top
    ft_singleplotER(cfg, cov);

    title(['sub-', sub, ': occipital channels'])
    set(gcf, 'color', 'w')

    cmap = cmocean('balance');

    f2 = figure;
    layout = '/project/2422120.01/BIDS/ActiCap_64Ch_DCC_customized.mat';
    cfg = [];
    cfg.figure = f2;
    cfg.interactive = 'no';
    cfg.layout = layout;
    cfg.parameter = 'avg';
    cfg.comment = 'no';
    cfg.colormap = cmap;
    cfg.colorbar = 'yes';
    cfg.zlim = 'maxabs';
    cfg.xlim = cov_time;
    ft_topoplotER(cfg, cov)

    sgtitle(sprintf('sub-%s: %.2f-%.2f sec', sub, cov_time(1), cov_time(2) ))
    set(gcf, 'color', 'w')

    %% Visualize cov matrix
%
%     f3 = figure;
%     imagesc(cov.cov);
%     c_lim = max(abs([min(cov.cov), max(cov.cov)]));
%     colorbar
%     caxis([-c_lim, c_lim]);
%     colormap(cmap)
%
    rank_cov = rank(cov.cov);
    rank_def = length(cov.label) - rank_cov;
%
%     title(['sub-', sub, ': rank deficient by ', num2str(rank_def)]);
%     set(gcf, 'color', 'w')

    %% beamform this

    cfg = [];
    cfg.method = 'lcmv';
    cfg.lcmv.fixedori = 'yes';
    cfg.sourcemodel = leadfield;
    cfg.elec = elec_aligned;
    cfg.lcmv.weightnorm = 'nai';
    cfg.headmodel = vol;
    cfg.lcmv.kappa = rank_cov;
    cfg.lcmv.keepfilter = 'yes';
    source = ft_sourceanalysis(cfg, cov);


    %% find all occipital voxels

    % boundaries in MNI coordinates: [x, y, z] by [min, max]
    bounds = [-84, 84;
        -120, -64;
        -12, 84] ./ 10;

    x_bool = bounds(1, 1) < template_grid.sourcemodel.pos(:, 1) & template_grid.sourcemodel.pos(:, 1) < bounds(1, 2);
    y_bool = bounds(2, 1) < template_grid.sourcemodel.pos(:, 2) & template_grid.sourcemodel.pos(:, 2) < bounds(2, 2);
    z_bool = bounds(3, 1) < template_grid.sourcemodel.pos(:, 3) & template_grid.sourcemodel.pos(:, 3) < bounds(3, 2);

    box_bool = [x_bool, y_bool, z_bool, template_grid.sourcemodel.inside];
    box_idx = find(all(box_bool, 2));

%% compute max voxel and plot time course

    [max_pow, max_idx] = max(source.avg.pow(box_idx));

    % find maximum timecourse and plot

    f5 = figure;
    plot(source.time, source.avg.mom{box_idx(max_idx)}); % initialize

    % add shaded area for cov time
    ylims = get(gca, 'ylim');
    ar = area(cov_time([1, 1, 2, 2]), [ylims, ylims([2, 1])], ...
        'handlevisibility', 'off');
    ar.FaceColor = [0.95, 0.95, 0.95];
    ar.EdgeColor = 'none';
    ar.BaseValue = ylims(1);
    hold on;

    % replot on top
    plot(source.time, source.avg.mom{box_idx(max_idx)});
    set(gca, 'xlim', [-0.3, 0.5]);
    title(['sub-', sub, ': time course in max. occip voxel'])
    set(gcf, 'color', 'w')

    %% Plot voxel box once

    % template MRI
    path_template_mri = fullfile(codetree.fieldtrip, ...
        '/template/anatomy/single_subj_T1.nii');

    % read template MRI (MNI space)
    temp_mri = ft_read_mri(path_template_mri);
    temp_mri.coordsys = 'mni';


    % plot the voxels to make sure:
     %if s_idx == 22 %??? ONLY PP 22?
     if sub == '022';
        source_demo = source;
        source_demo.pos = template_grid.sourcemodel.pos;
        source_demo.avg.pow(:) = 0;
        source_demo.avg.pow(box_idx) = 1;


        cfg = [];
        cfg.parameter = 'avg.pow';
        cfg.interpmethod = 'nearest';
        source_demo_int = ft_sourceinterpolate(cfg, source_demo, temp_mri);

        fx = figure;
        cmap = colormap(cmocean('ice'));

        % plot voxel box
        cfg = [];
        cfg.figure = fx;
        cfg.funcolormap = cmap;
        cfg.method = 'ortho';
        cfg.funparameter = 'pow';
        cfg.funcolorlim = 'zeromax';
        ft_sourceplot(cfg, source_demo_int);

    end

    %% Interpolate and plot

    % atlas
    path_aal_atlas = fullfile(codetree.fieldtrip, ...
        'template/atlas/aal/ROI_MNI_V4.nii');
    atlas = ft_read_atlas(path_aal_atlas);


    source.pos = template_grid.sourcemodel.pos;

    % interpolate power on template brain
    cfg = [];
    cfg.parameter = 'avg.pow';
    cfg.interpmethod = 'nearest';
    source_int = ft_sourceinterpolate(cfg, source, temp_mri);

    % utilities
    cmap = colormap(cmocean('thermal'));

    % plot power
    f4 = figure;
    cfg = [];
    cfg.figure = f4;
    cfg.atlas = atlas;
    cfg.funcolormap = cmap;
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0, max_pow];  % saturate occipitally
    cfg.location = source.pos(box_idx(max_idx), :) .* 10; % convert to mm
    ft_sourceplot(cfg, source_int);

    title(sprintf('sub-%s: %.2f-%.2f sec', sub, cov_time(1), cov_time(2) ))
    set(gcf, 'color', 'w')


%% arrange figures (spreadfigures need to be present on path)

    %spreadfigures([f1, f2, f4, f5], 'horizontal');
%Spreadfigures([f1, f2, f4, f5], 'square');


%end % enable if running multiple pp
