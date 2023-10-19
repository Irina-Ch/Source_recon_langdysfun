%% MRI Preprocessing for EEG Forward modelling

% Possible preprocessing steps for MRIs prior to EEG forward modelling
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script does optional preprocessing on the MRI before segmenting the
% MRI for BEM modelling.

%% Set option

% Optional bias field correction of the MRI
bias_field_corr = 0;   

%% Attention!

% This pipeline assumes nifti files. If your MRI is in DICOM format, you
% can convert it to nifti by reading it in with FieldTrip and
% then writing it out again using ft_write_mri() - see below for an example

% IMPORTANT! SPM cannot deal with zipped nifti files - i.e., the ending
% should be .nii and not .nii.gz

% IMPORTANT! Do not use any other programs for converting to nifti, unless
% you are very sure they give you good output - some programs might by
% default strip the skull or do other processing which will pose problems
% for this pipeline.

%% Example of converting to nifti:
try
  mri_tmp = ft_read_mri(projpath.mri); % if there is no nifty in the folder, save as one:
catch
  disp('There is no nifti available!') 
  n = dir([projpath.dicoms, '/*.IMA']); % get the names in the directory
    mri_tmp = ft_read_mri([projpath.dicoms, '/', n(1).name]); % n(1).name contains first .ima name
    ft_write_mri(projpath.mri, mri_tmp, 'dataformat', 'nifti');
  clear mri_tmp
end
% for a series of DICOM files, provide the name of any of the files in the series (e.g. the first one). 
% The files corresponding to the whole volume will be found automatically.

%mri_tmp = ft_read_mri('/project/2422120.01/BIDS/sub-016/anat/dicom_s016/00001_1.3.12.2.1107.5.2.19.45416.2023052311482376728700069.IMA')
%mri_tmp = dicomread([projpath.dicoms, '/00001_1.3.12.2.1107.5.2.19.45416.2023070313285686598602306.IMA']);
%dicom_info = dicominfo([projpath.dicoms, '/00001_1.3.12.2.1107.5.2.19.45416.2023070313285686598602306.IMA']);

% Optionally, you can also read the DICOM metadata

% projpath.mri is: '/project/2422120.01/BIDS/sub-001/anat/sub-001_T1w.nii'

% NOW MAKE SURE to go and update projpath.mri to point to the new nifti 
% in ldf_00_setup.m and run ldf_00_setup.m again before continuing!

%% Plot the MRI 

mri_in = ft_read_mri(projpath.mri);
ft_sourceplot([], mri_in);

%% Bias field correction

% MRI scans can be corrupted by a bias (introduced by the MRI machine)
% which makes segmentation algorithms fails. The bias field signal basically 
% changes the grey-value of voxels across space, which poses difficulties
% for the segmentation algorithms that rely on the grey value / contrast to
% classify voxels into tissue types.
% Here we correct for this by running a bias field correction in SPM. 
% This should be done if the segmentation fails (visible holes in the scalp 
% layer) or if the scalp mesh does not close.

% This wil take some time, so be prepared to wait!

if bias_field_corr

    addpath(toolboxes.spm12);  %#ok
    
    % set the configurations
    % some of these configuration values are taken over from 
    % https://layerfmri.com/2017/12/21/bias-field-correction/
    spm_config{1}.spm.spatial.preproc.channel.vols = {[projpath.mri, ',1']};
    spm_config{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    % a FWHM of 30 is the lightest option SPM supports. Change this to 40
    % if you still experience problems with segmentation:
    spm_config{1}.spm.spatial.preproc.channel.biasfwhm = 30; % 30; NB! this one, try 40, go up more
    spm_config{1}.spm.spatial.preproc.channel.write = [1 1];
    spm_config{1}.spm.spatial.preproc.warp.affreg = 'mni';
    spm_config{1}.spm.spatial.preproc.warp.cleanup = 1;
    spm_config{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    spm_config{1}.spm.spatial.preproc.warp.fwhm = 0;
    spm_config{1}.spm.spatial.preproc.warp.samp = 3;
    spm_config{1}.spm.spatial.preproc.warp.mrf = 1;
    spm_config{1}.spm.spatial.preproc.warp.write = [0 0];
    
    % run the job - it will automatically save the bias field corrected MRI
    % (and some by-products) to the MRI folder.
    spm('defaults', 'FMRI')
    spm_jobman('initcfg');
    spm_jobman('run', spm_config);

    % if we do this, we have to reload the bias-field-corrected MRI now
    % since we want to continue with this.
    mri_in = ft_read_mri(projpath.mri_bfc);
    ft_sourceplot([], mri_in);  % check that it looks decent

end

%% reslice MRI

% Reslicing an MRI makes sure that the voxels are isotropic - i.e., the
% edges of the voxels have the same length in all directions.

% This also saves the MRI, which is crucial for further processing.
%NB! Do NOT RESCICE FOR SEGMENTATION PROBLEMATIC PPS

%% Get non-resliced MRI for problematic pps

if str2num(id) ~= 2 & str2num(id) ~= 17 & str2num(id) ~= 6 & str2num(id) ~= 7 & str2num(id) ~= 8 & ...
        str2num(id) ~= 9 & str2num(id) ~= 11 & str2num(id) ~= 12 % for problematic segmentation pps
    mri_resl = ft_volumereslice([], mri_in); % Ir: TYPO changed to 'mri_in' from 'mri'
    save(projpath.mri_resl, 'mri_resl'); 
end

%%% from some of V's scripts: WE ALREADY RESLICE, though
% convert MRI to ctf coordinates and make sure resolution is mm (reslicing):
%     cfg                     = [];
%     cfg.method              = 'interactive';
%     cfg.coordsys            = 'ctf'; 
%     mri_aligned             = ft_volumerealign(cfg,mri);
%     cfgres.resolution       = 1; % in mm
%     mri_reslice             = ft_volumereslice(cfgres, mri_aligned); 
