%% Align MRI to defined space

% Align MRI roughly to BTI space for successful segmentation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the right MRI

% It is important that you now choose the preprocessed MRI that you want to use
% for creating your head model (i.e. the resliced and potentially bias 
% corrected one). 
% Otherwise, there will be a mismatch between the coordinate 
% systems and your electrode positions will not align with the head model.

% For the following pipeline to be successful, reslicing your MRI is highly
% recommended. Bias field correction should be chosen in case it is
% necessary (i.e. if the segmentation failed).

if str2num(id) == 2 | str2num(id) == 17 | str2num(id) == 6 | str2num(id) == 7 | str2num(id) == 8 | ...
        str2num(id) == 9 | str2num(id) == 11 | str2num(id) == 12 % for problematic segmentation pps
    mri_path = fullfile(projpath.base, mri_name); % to use non-resliced mri 
else
    mri_path = projpath.mri_resl;
end

% find and load the right MRI
if strcmp(mri_path, projpath.mri)
    warning('You should probably reslice and realign your MRI!')
    mri = ft_read_mri(mri_path);
elseif strcmp(mri_path, projpath.mri_resl)
    load(projpath.mri_resl)
    mri = mri_resl; clear mri_resl
elseif strcmp(mri_path, projpath.mri_bfc_resl)
    load(projpath.mri_resl)
    mri = mri_resl; clear mri_resl
else
    error('Do not know how to load specified MRI');
end

%% Identify the fiducials on the MRI - NOT COREGISTRATION STEP!!!

% Now we want to identify the fiducial points on the MRI. These points
% we'll use to align the MRI and the EEG electrodes. So be careful here -
% this directly influences the accuracy of your forward and thus your
% inverse model!

% Follow the instructions that get printed to the command window and mark:
% Nasion - n , LPA - l, RPA - r, and a z-point - z.
% You do not have to be *overly* precise here - this is NOT THE
% COREGISTRATION!
% We just bring the MRI into CTF space which will help us with segmentation
% as well as with using the structural scan for the coregistration step.

cfg = [];
cfg.method   = 'interactive';
cfg.coordsys = 'bti';
mri_aligned = ft_volumerealign(cfg, mri);

%% Save output

% we save the fiducials in case we ever need them again (e.g. if we have to
% coregister again)
fids = mri_aligned.cfg.fiducial;
fids.coordsys = mri_aligned.coordsys;
fids.mri2meg = mri_aligned.transform;
fids.origtfm = mri_aligned.transformorig;

% this is a manual step, so we save the result!
save(projpath.fids, 'fids');

% And more importantly, we save the aligned MRI, too.
% From now on, we will only work with this aligned MRI!
save(projpath.mri_aligned, 'mri_aligned');
