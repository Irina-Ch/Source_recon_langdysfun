function [source_difference] = getSOURCEnoun (sbjnumb)

% Takes raw data, cleans it, redefines to picture onset and produces source
% estimation calcualtions for NOUN task. Use this function for parallel
% processing.

%% Get all file names
task = 2; % NOUN TASK

id = sbjnumb;

projpath = []; 
projpath.base = ['/project/2422120.01/BIDS/sub-', id, '/anat'];
projpath.vol = fullfile(projpath.base, sprintf('%s_vol.mat', id));
projpath.fwd = fullfile(projpath.base, sprintf('%s_fwd.mat', id));
projpath.elec_aligned = fullfile(projpath.base, sprintf('%s_elec_aligned', id));


%% Load files
load(projpath.fwd) % the forward model
load(projpath.vol) % the headmodel again (to plug into ft_sourceanalysis)
load(projpath.elec_aligned) % load aligned electrodes: use elec_realigned

% get clean EEG data:
sbjnumb = id(2:3); %two last digits from id
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
%cfg.foilim          = [8 30]; 
cfg.foi = 10;
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
cfg              = [];
cfg.method       = 'dics';
cfg.sourcemodel  = leadfield;
cfg.headmodel    = vol;
cfg.elec = elec_realigned; 
cfg.frequency    = 10; % single number
cfg.dics.fixedori = 'yes';
cfg.dics.weightnorm = 'unitnoisegain'; % normalise here (since no control/baseline) 
%cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';  % We want to reuse the calculated filter later on
source_all = ft_sourceanalysis(cfg, powcsd_all);

%% Apply common filter to both conditions
%run ft_sourceanalysis again - for both conditions - but specify that we want to use the filter we just computed.
cfg              = [];
cfg.method       = 'dics';
cfg.elec = elec_realigned; 
cfg.sourcemodel  = leadfield;
cfg.headmodel    = vol;
%cfg.frequency    = 10; 
cfg.sourcemodel.filter = source_all.avg.filter;  % We apply the previously computed spatial filter

source_con  = ft_sourceanalysis(cfg, powcsd_con);
source_unc = ft_sourceanalysis(cfg, powcsd_unc);


%% Compute the difference between conditions
source_difference = source_con;
source_difference.avg.pow = (source_con.avg.pow - source_unc.avg.pow) ./ ...
                      ((source_con.avg.pow + source_unc.avg.pow)/2);


