function [mri_fids, surface_fids] = identify_points(mri_mesh, scan_mesh)
%
% Function to identify identical face points on two meshes. Plots the two
% meshes and prompts the user to identify points on them. Returns those
% points from both meshes.
%
% Parameters:
% mri_mesh: the mesh from the MRI (scalp)
% scan_mesh: the mesh from the structural scan
%
% AUTHOR: Britta U. Westner <britta.wstnr@gmail.com>
% CREDIT: Based on code from Sarang S. Dalal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure init
fig_1 = figure(1);

% plot structural scan
ax_scan = subplot(1, 2, 1);
trisurf(scan_mesh.faces, ...
    scan_mesh.vertices(:, 1), ...
    scan_mesh.vertices(:, 2), ...
    scan_mesh.vertices(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'red');
ax_scan.Interactions = dataTipInteraction;
disableDefaultInteractivity(ax_scan);  % turn off irritating data tips
axis equal
view([-90, -25])  
%lightangle(-125, -50); % ORIGINAL; change here per pp
%lightangle(0, -170); % sub022
lightangle(0, -190); % sub002, sub005; sub013; sub014; sub024; sub025

% plot MRI
ax_mri = subplot(1,2,2);

trisurf(mri_mesh.faces, ...
    mri_mesh.vertices(:, 1), ...
    mri_mesh.vertices(:, 2), ...
    mri_mesh.vertices(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'white');
ax_mri.Interactions = dataTipInteraction;
disableDefaultInteractivity(ax_mri);
axis equal
view([90, 0]);
lightangle(45,-50);

fid_prompts = {'Please select middle of left eye.', ...
    'Please select middle of right eye.', ...
    'Please select nasion (bridge of nose).', ...
    'Please select tip of nose.', ...
    'Please select left corner of nose.', ...
    'Please select right corner of nose.'};

% iterate over the points to identify
surface_fids = zeros(length(fid_prompts), 3);
mri_fids = surface_fids;
for ii=1:length(fid_prompts)
    
    % structural scan
    axes(ax_scan);  %#ok
    h_title = title(fid_prompts{ii});
    pause
    surface_fids(ii,:)=select3d;
    delete(h_title);

    % MRI 
    axes(ax_mri)  %#ok
    h_title = title(fid_prompts{ii});
    pause
    mri_fids(ii,:)=select3d;
    delete(h_title)
    
end

pause(1);  % make closing less promptly
close(fig_1);
