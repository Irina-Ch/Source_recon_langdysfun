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
% mri_tmp = ft_read_mri(projpath.mri);
% ft_write_mri('~/MEG/test_mri/001-T1.nii', mri_tmp, 'dataformat', 'nifti')

% NOW MAKE SURE to go and update projpath.mri to point to the new nifti 
% in ldf_00_setup.m and run ldf_00_setup.m again before continuing!

%% Plot the MRI 

% This will also make sure the file actually exists :)

mri_in = ft_read_mri(projpath.mri);
ft_sourceplot([], mri_in);

%% Bias field correction

% MRI scans can be corrupted by a bias (introduced by the MRI machine)
% which makes segmentation algorithms fails. The bias field signale basically 
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
    spm_config{1}.spm.spatial.preproc.channel.biasfwhm = 30;  
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
mri_resl = ft_volumereslice([], mri_in); % changed to 'mri_in' from 'mri'
save(projpath.mri_resl, 'mri_resl'); 


