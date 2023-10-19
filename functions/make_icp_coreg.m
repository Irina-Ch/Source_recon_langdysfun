function [tfm, scan_mesh] = make_icp_coreg(mri_mesh, scan_mesh, face_points, tfm_0)
%
% Function that uses pre-selected face points and reduces the meshes to
% the face only. The face meshes are then used to coregister the structural
% scan and the MRI via ICP. Returns the transform from structural scan to
% MRI and the structural scan mesh including the updated vertices.
%
% Parameters:
% mri_mesh: the mesh from the MRI (scalp)
% scan_mesh: the mesh from the structural scan
% face_points: 4 coordinate points that identify the boundaries of the face
% tmf_0: the transform obtained from a first, coarse coregistration
%
% AUTHOR: Britta U. Westner <britta.wstnr@gmail.com>
% CREDIT: Based on code from Sarang S. Dalal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use the face points to restrict the mesh to the face

min_x = mean([face_points(1, 1); face_points(4, 1)]);
min_y = face_points(1, 2);
max_y = face_points(4, 2);
min_z = face_points(3, 3);
max_z = face_points(2, 3);

% find the vertices to keep
del_x = find(scan_mesh.vertices_0(:, 1) < min_x);
del_y = find(scan_mesh.vertices_0(:, 2) < min_y | scan_mesh.vertices_0(:, 2) > max_y);
del_z = find(scan_mesh.vertices_0(:, 3) < min_z | scan_mesh.vertices_0(:, 3) > max_z);

del_vert_faces = unique([del_x; del_y; del_z]);

del_x = find(mri_mesh.vertices(:, 1) < min_x);
del_y = find(mri_mesh.vertices(:, 2) < min_y | mri_mesh.vertices(:, 2) > max_y);
del_z = find(mri_mesh.vertices(:, 3) < min_z | mri_mesh.vertices(:, 3) > max_z);

del_vert_mri = unique([del_x; del_y; del_z]);

scan_face = scan_mesh;
scan_face.vertices_0(del_vert_faces) = NaN;

mri_face = mri_mesh;
mri_face.vertices(del_vert_mri) = NaN;

%% Plot to double check

fig1 = figure(1);
set(fig1, 'toolbar', 'none')
subplot(1, 2, 1);
trisurf(scan_face.faces, ...
    scan_face.vertices_0(:, 1), ...
    scan_face.vertices_0(:, 2), ...
    scan_face.vertices_0(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'red');
axis equal
view([90, 0]); 
lightangle(-125, -50);
title('Struct. scan')

subplot(1,2,2);
trisurf(mri_face.faces, ...
    mri_face.vertices(:, 1), ...
    mri_face.vertices(:, 2), ...
    mri_face.vertices(:, 3), ...
    'EdgeColor', 'none', 'FaceColor', 'white');
axis equal
view([90, 0]);
lightangle(-125, -50);
title('MRI');

subtitle('These are the 3D data used for ICP'); %Ir orig: suptitle gives error

%% Reduce coreg error via ICP of facial points

sel_mri = setdiff(1:length(mri_mesh.vertices), del_vert_mri);
sel_scan = setdiff(1:length(scan_mesh.vertices_0), del_vert_faces);

% ICP
[tr, tt] = icp(unique(mri_mesh.vertices(sel_mri,:), 'rows')', ...
    unique(scan_mesh.vertices_0(sel_scan,:), 'rows')');

% compute transformation matrix and apply it
tfm = [tr, tt; 0 0 0 1] * tfm_0;
scan_mesh.vertices_1 = nut_coordtfm(scan_mesh.vertices, tfm);

