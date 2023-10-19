%% Identifying the electrodes on a structural scan

% Identify the fiducials on the MRI for coregistration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Britta U. Westner <britta.wstnr[at]gmail.com>
% LICENCE: GNU General Public License v3.0
% script is based on:
% https://www.fieldtriptoolbox.org/tutorial/electrode/
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the procedure to follow in case you work with a structural scan
% of your participant and the EEG cap.
% It is strongly recommended to also read:
% https://www.fieldtriptoolbox.org/tutorial/electrode/

%% Set paths

% Since this might be a step not everyone follows, we set the path to the
% necessary ingredients here. If you will use this script, you might
% consider moving this to your setup file (your ldf_00_setup.m or wherever you 
% store the paths in your project).

% this is the path to your structural scan:
projpath.struct_scan = fullfile(projpath.base, sprintf('%s_model.obj', id));

% path to the template for the EEG cap layout
% IMPORTANT: if you work at the DCC, then there is a specific layout 
% "ActiCap_64Ch_DCC_customized.mat" that you can find in our common
% folders.
template_grid = load(fullfile(toolboxes.fieldtrip, ...
    '/template/layout/acticap-64ch-standard2.mat'));

%% Load the data and convert units

% we read the structural scan as a "headshape"
head_surface = ft_read_headshape(projpath.struct_scan);

% make sure things are in millimeters:
head_surface = ft_convert_units(head_surface, 'mm'); 

%% Look at the model surface

figure;
ft_plot_mesh(head_surface); 

%% Identify electrode locations

% Now we identify the electrodes on the structural scan. We identify them
% on the original structural scan mesh - that is the one we also use for
% the coregistration later. You can skip identifying the fiducials. 
% Here it is important to match the electrodes identified and the electrode
% names correctly (otherwise you will end up with electrodes that have the
% coordinates of other electrodes!).
% The function get_electrode_order() supplies you with an order that you
% can follow in the identification process. You can also create your own
% order - but you have to be careful that you supply this order as channel
% names as well then!

% This function to get the electrode order is shared in the subfolder
% /functions
elec_order = get_electrode_order('acticap_dcc'); %NB! I also have lpa,nas, rpa marked

% Important: if FieldTrip asks you if you want to change the anatomical
% labels for the axes, answer NO (n).

cfg = [];
cfg.method = 'headshape';
cfg.channel = elec_order;
elec = ft_electrodeplacement(cfg, head_surface);

%% Save File

save(projpath.elec, 'elec');

