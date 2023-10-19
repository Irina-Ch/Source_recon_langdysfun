%% Check leadfields CTXTVerb
%addpath('/home/langdysfun/irichu/Documents'); % Spreadfigures.m
toolboxes.fieldtrip = '/home/langdysfun/irichu/Documents/toolboxes/fieldtrip/';

% Fieldtrip needs to be on path.
codetree.cmocean = fullfile('/home/langdysfun/irichu/Documents/cmocean-main');
addpath(codetree.cmocean)

%% plot the leadfield

pos = [18, -98, 6] / 10; % convert to cm

% load the template
template_grid = load(fullfile(toolboxes.fieldtrip, ...
    '/template/sourcemodel/standard_sourcemodel3d10mm.mat'));

% find the position index
pos_idx = dsearchn(template_grid.sourcemodel.pos, pos);

%% Run across subjects

%subjects = arrayfun(@(X) sprintf('%.3d', X), 1:28, 'uniformoutput', 0)';
%subjects = setdiff(subjects, {'002', '006', '007', '008', '009', '011','012', '017'}); %remove these
subjects = {'002', '006', '007', '008', '011', '012', '017'};
% s009 still not working

for s_idx = 1:length(subjects)

    sub = subjects{s_idx};

    bids_path = ['/project/2422120.01/BIDS/sub-', sub, '/anat'];
    load(fullfile(bids_path, [sub, '_fwd.mat']));

    %% Visualize the position once

    if s_idx == 1
        inside_idx = find(leadfield.inside);

        figure;
        plot3(leadfield.pos(inside_idx, 1), leadfield.pos(inside_idx, 2), ...
            leadfield.pos(inside_idx, 3), '.', 'color', [0.73, 0.73, 0.73]); hold on
        plot3(leadfield.pos(pos_idx, 1), leadfield.pos(pos_idx, 2), ...
            leadfield.pos(pos_idx, 3), '.', 'color', [1, 0, 0], 'markersize', 15)
        axis equal
        set(gcf, 'color', 'w')
    end

    %% plot this

    figure;
    for ori = 1:3

        subplot(1, 3, ori)
        data.label = leadfield.label;
        data.avg = leadfield.leadfield{pos_idx}(:, ori);
        data.time = 1;

        cmap = cmocean('curl');

        layout = '/project/2422120.01/BIDS/ActiCap_64Ch_DCC_customized.mat';
        cfg = [];
        cfg.zlim = [-1.2, 1.2]*10^-4;
        cfg.figure = 'gca';
        cfg.interactive = 'no';
        cfg.layout = layout;
        cfg.parameter = 'avg';
        cfg.comment = 'no';
        cfg.colormap = cmap;
        ft_topoplotER(cfg, data)

    end

    set(gcf, 'color', 'w')
    sgtitle(['Leadfields: ', sub])
end

% manage plots; needs speadfigures on path:
Spreadfigures([figs{:}], 'square')