%% Testing solutions for segmentation issues

%addpath('/home/langdysfun/irichu/Documents'); % Spreadfigures.m
dirfig = ['/project/2422120.01/BIDS/sub-', id, '/anat/figures/']; % for saving figures


test_mri = ft_read_mri(projpath.mri); % run ldf_00
test_mri.coordsys = 'spm';

%%

cfg = [];
cfg.output = {'brain', 'skull', 'scalp'};
cfg.spmversion = 'spm12';
cfg.spmmethod = 'old';
seg = ft_volumesegment(cfg, test_mri);

save(projpath.seg, 'seg')
%%
load(projpath.seg)

%% Use private FieldTrip function to attempt filling holes
% we need to jump directories because you can't put a private function
% directly on your path in MATLAB !

seg_fill = seg;  % make a copy of seg

% jump directory as dir with function is private and thus cannot be added to path
pw_dir = pwd;
[~, ft_path] = ft_version;
cd(fullfile(ft_path, 'private'));

% fill along z axis
seg_fill.scalp = volumefillholes(seg.scalp, 3);

% for sub006:
% seg_fill.skull = volumefillholes(seg.skull, 3);


% jump back to previous directory:
cd(pw_dir);

%% Dilate the scalp segmentation

seg_dil = seg_fill;  % work with filled segmentation

% dilation:
se_scalp = strel('cube', 3); % s009: 2 pixels do not work; strel("cube",w) creates a 3-D cubic structuring element whose width is w pixels.
% strel(): A strel object represents a flat morphological structuring element, which is an essential part of morphological dilation and erosion operations.
seg_dil.scalp = imdilate(seg_fill.scalp, se_scalp); % grayscale or binary dilation of an image
                           %image          % structuring element

%% SKIP this step for s006 and s009 !
% correct dilation by subtracting the skull:
seg_dil.scalp = seg_dil.scalp - seg_dil.skull;
% UPDATE: this might be better left to the meshing step and not done here.
% At least for MRI006 and MRI009, this seems to end up creating holes again.
% Note that not doing this step could also break things ... So if you omit
% it and run into trouble with the headmodel, it might be worth getting it
% back into the mix.
% seg_dil.scalp(seg_dil.skull) = 0;

%% Also add de-islanding for worst cases: s006
% for these to run, iso2mesh should be on path (cfg.method = 'iso2mesh';)
% e.g. from ft_prepare_mesh

if(strcmp(id, '006'))
    % for MRI 006, but not MRI 009

    new_scalp = deislands3d(seg_dil.scalp);
    seg_dil.scalp(new_scalp) = 1; % correct scalp using the deislanding
%     seg_dil.scalp(seg_dil.skull) = 0; % correct with skull

    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'scalp';
    cfg.colorbar = 'no';  % no functional data
    ft_sourceplot(cfg, seg_dil);
    sgtitle(['sub-',id, ' after deislanding'])
    temp = [dirfig 'after_deislanding_sub-' id '.png']; print('-dpng',temp);
end

%% Go on, mesh both and plot to compare
% in the end we get 'mesh' based on dilated segmentation: save it!
segments = {seg, seg_dil};
seg_names = {'original segm.', 'dilated segm.'};
figs = cell(6, 1); %save figures in grid
count = 0;

for jj = 1:length(segments)

    cfg = [];
    cfg.tissue = {'brain', 'skull', 'scalp'};
    cfg.method = 'iso2mesh';
    cfg.spmversion = 'spm12';
    cfg.numvertices = 1e4;
    mesh = ft_prepare_mesh(cfg, segments{jj});

    % Downsample and repair the meshes

    target_sizes = [2000, 2000, 2000];

    for ii = 1:length(mesh)
        [mesh(ii).pos, mesh(ii).tri] = meshresample(mesh(ii).pos, ...
            mesh(ii).tri, ...
            target_sizes(ii)/...
            size(mesh(1).pos, 1));
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
            mesh(ii).tri, 'dup');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
            mesh(ii).tri, 'isolated');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
            mesh(ii).tri, 'deep');
        [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos, ...
            mesh(ii).tri, 'meshfix');
    end

    % Plot the meshes
    views = [0, 90, 0; 90, 0, 0; -90, 0, 0];  % three views

    for ff = 1:3

        count = count + 1;
        figs{count} = figure;
        ft_plot_mesh(mesh(1), 'facecolor', 'red');
        hold on; % pause
        ft_plot_mesh(mesh(2), 'facealpha', 0.25, 'facecolor', 'blue');
        hold on; % pause
        ft_plot_mesh(mesh(3), 'facealpha', 0.25, 'edgecolor', [0.8, 0.8, 0.8]);
        title(sprintf('%s: %s', 'sub-',id, seg_names{jj}))

        view(views(ff, :))
    end
%    % attempt the headmodel:
%     try
%         cfg = [];
%         cfg.method = 'dipoli';  % dipoli, bemcp, or openmeeg (the latter needs installation)
%         cfg.conductivity = [0.33 0.0041 0.33];
%         vol = ft_prepare_headmodel(cfg, mesh);
%     end

end

% manage plots; needs speadfigures on path:
Spreadfigures([figs{:}], 'square')

%plot the dialted mesh:
    figure;
    ft_plot_mesh(mesh(1), 'facecolor', 'red');
    hold on; pause
    ft_plot_mesh(mesh(2), 'facealpha', 0.25, 'facecolor', 'blue');
    hold on; pause          
    ft_plot_mesh(mesh(3), 'facealpha', 0.25, 'edgecolor', [0.8, 0.8, 0.8]);
    view([0 90 0])
    
    temp = [dirfig 'mesh_sub-' id '.png']; print('-dpng',temp);
    

% save the mesh:
save(projpath.mesh, 'mesh', '-v7.3');

