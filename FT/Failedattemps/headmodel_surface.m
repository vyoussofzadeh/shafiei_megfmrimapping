
%% construct mesh from brain
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;
mesh_brain = ft_prepare_mesh(cfg, individual_seg);

%%
% MEG
cfg = [];
cfg.method = 'singleshell';
headmodel_mne_meg = ft_prepare_headmodel(cfg, mesh_brain);

%%
figure; hold on
% ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_vol(headmodel_mne_meg, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;

ft_plot_mesh(surface_grid, 'edgecolor', 'k'); camlight

%%
% Get transformation matrix from voxel to neuromag-coordsys
transform_vox2neuromag = individual_seg.transform;
transform_vox2neuromag = mri_realigned.transform;
% Get transformation matrix from spm to neuromag
T = transform_vox2neuromag/transform_vox2spm_rs;

sourcemodelT = ft_transform_geometry(mri_realigned.transform, surface_grid);

%%
headmodel_mne_meg = ft_convert_units(headmodel_mne_meg, 'mm');
sourcemodelT = ft_convert_units(sourcemodelT, 'mm');

%%
figure; hold on
ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight

%%
figure; hold on
ft_plot_vol(individual_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(individual_grid, 'edgecolor', 'k'); camlight

%%
% cfg          = [];
% % cfg.method   = 'interactive';
% cfg.method = 'fiducial'; % the following voxel coords were determined interactive
% cfg.coordsys = 'acpc';
% cfg.fiducial.ac   = fid.NCS.AC;
% cfg.fiducial.pc   = fid.NCS.PC;
% cfg.fiducial.xzpoint  = fid.NCS.IH;
% cfg.fiducial.right   = fid.SCS.RPA;
% mri_spm     = ft_volumerealign(cfg, mri_realigned);
% 
% cfg = [];
% cfg.method = 'fiducial';
% cfg.coordsys = 'neuromag';
% cfg.fiducial.nas   = fid.SCS.NAS;
% cfg.fiducial.lpa   = fid.SCS.LPA;
% cfg.fiducial.rpa   = fid.SCS.RPA;
% cfg.spmversion     = 'spm12';
% mri_neuromag     = ft_volumerealign(cfg, mri_realigned);

% 
% transform_vox2acpc = mri_acpc.transform;
% transform_vox2neuromag = mri_realigned.transform;
% transform_acpc2neuromag = transform_vox2neuromag/transform_vox2acpc;
% 
% 
% %%
% sourcemodel = ft_transform_geometry(transform_acpc2neuromag, surface_grid);
% individual_headmodel1 = ft_transform_geometry(transform_acpc2neuromag, individual_headmodel);

%%
% figure; hold on
% ft_plot_vol(individual_headmodel1, 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);

%%
% surface_sourcemodel = ft_convert_units(surface_sourcemodel, 'mm');
% figure; ft_plot_mesh(surface_sourcemodel, 'edgecolor', 'k'); camlight
% hold on
% ft_plot_vol(headmodel_mne_meg, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% ft_plot_mesh(individual_grid, 'edgecolor', 'k'); camlight

%%
% from voxel to spm-coordsys.
% load transform_vox2spm.mat
% transform_vox2spm_rs = mri_spm.transform
% % Load segmented mri and get transform
% 
% % Get transformation matrix from voxel to neuromag-coordsys
% transform_vox2neuromag = mri_neuromag.transform
% 
% % Get transformation matrix from spm to neuromag
% T = transform_vox2neuromag/transform_vox2spm_rs
% % T = mri_realigned.transformorig/mri_realigned.transform;
% T = mri_realigned.transform;

%%
% sourcemodelT = ft_transform_geometry(transform_vox2spm_rs, surface_sourcemodel);
% sourcemodelT2 = ft_transform_geometry(transform_vox2neuromag, sourcemodelT1);

%%
transform = [
    0 -1 0 0;
    1 0 0 0;
    0 0 1 0;
    0 0 0 1
    ];
%%
sourcemodelT = ft_transform_geometry(transform, surface_sourcemodel);

%%
% cfg = [];
% % cfg.method    = 'interactive';
% cfg.method = 'headshape';
% cfg.headshape.interactive = 'yes';
% cfg.headshape.icp = 'yes';
% cfg.coordsys = 'neuromag';
% cfg.headshape      = surface_sourcemodel;
% % cfg.headshape = headmodel_mne_meg.bnd(1);
% mri_realignedT = ft_volumerealign(cfg,mri_realigned);

%%
figure; hold on
ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight
view([90,0])
view([-90,0])
view([0,0])

