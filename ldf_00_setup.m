%% Setup file

% Setup for forward modelling
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this file, we set the paths we will need later and load toolboxes we
% will use. Make sure to update the paths!

%% Paths for toolboxes

% Here we just store the toolbox paths for easy access - change the paths
% to wherever you have saved those toolboxes. If you do not have them
% installed yet, you can find info on where to obtain them at the very end
% of this script!

toolboxes.fieldtrip = '~/Documents/code_dev/fieldtrip';
toolboxes.spm12 = '~/Documents/code_dev/spm12';  % needed for bias field corr.
toolboxes.nifti = '~/Documents/toolboxes/NIfTI_20140122/';  % needed for reslicing

%% Paths for data 

projpath = [];  % initialize

% first, the path to the folder of the MRI, change this accordingly:
% SPM is particular about the way the path is set, so if you plan to use
% SPM - e.g. the bias field correction step - make sure to not use the
% tilde shortcut for doing so (i.e. '~/MEG/test_mri' below), but give the
% full path here
projpath.base = '/home/predatt/briwes/MEG/test_mri/';

% second, the name of the MRI and a participant ID, also change:
% This should be the path to a nifti file - if your data is DICOM or else,
% check ldf_01_preprocess_mri.m for some info on how to convert.
mri_name = '002_T1w.nii';
id = '002';

% The rest of this section is automatic - no need to change further paths.
% But make sure to read the comments to know what's gonna happen and where
% data will be saved.

% generate full path to MRI out of the information
projpath.mri = fullfile(projpath.base, mri_name);  % full path to MRI

% We generate the other MRI file names, don't change this, as it partly
% relies on how SPM saves the data
projpath.mri_bfc = fullfile(projpath.base, ['m', mri_name]);  
mri_stem = split(mri_name, '.');
projpath.mri_resl = fullfile(projpath.base, [mri_stem{1}, '_resl.', mri_stem{2}]);  
projpath.mri_bfc_resl = fullfile(projpath.base, ['m', mri_stem{1}, '_resl.', mri_stem{2}]);  

% Now we generate some filenames we can use for saving and reusing
% generated output.
% This will assume that you want to save all the data that is being generated
% into the same folder. That should be fine for this demo, but maybe not be
% the smartest choice when you do extensive data analysis.
projpath.seg = fullfile(projpath.base, sprintf('%s_seg.mat', id));
projpath.vol = fullfile(projpath.base, sprintf('%s_vol.at', id));
projpath.grid = fullfile(projpath.base, sprintf('%s_grid.at', id));
projpath.fwd = fullfile(projpath.base, sprintf('%s_fwd.mat', id));


%% Add FieldTrip to path

% This is my way of adding FieldTrip to my path (I also like to check for the 
% Github branch I am on as a safety precaution). You can just replace this
% whole section by however you add FieldTrip to your path usually, e.g. by 
% simply doing:
% addpath(toolboxes.fieldtrip)

git.branch = 'master';
git.head = ['ref:refs/heads/', git.branch];

try
    ft_defaults
catch
    warning('Fieldtrip is not on your path yet, adding it.');
    read_git = fscanf(fopen(fullfile(toolboxes.fieldtrip, '.git/HEAD'), 'r'), '%s');

    if strcmp(git.head, read_git) == 0
        error('Wrong git HEAD. Change git branch to %s.', git.branch)
    end
    addpath(toolboxes.fieldtrip)
    ft_defaults
end

[ft_ver, ft_path] = ft_version;
fprintf('You are using Fieldtrip on path %s, branch %s, version hash %s \n', ...
                ft_path, git.branch, ft_ver);

%% Information on necessary toolboxes

% FieldTrip
% Download information here:
% https://www.fieldtriptoolbox.org/download/

% SPM12
% SPM version 12 can be obtained here:
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% make sure to choose version 12!

% reslice_nii as part of Nifti:
% https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
