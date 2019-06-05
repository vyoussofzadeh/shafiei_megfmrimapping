% template_mri = ft_read_mri('E:\My Matlab\SPM\spm12_4\spm12\canonical\avg152T1.nii'); % based on mni151
template_mri = ft_read_mri('F:\My Matlab\SPM\spm12_4\spm12\spm12\canonical\avg152T1.nii'); % based on mni151
template_mri.coordsys = 'spm';  % inform fieldtrip that 'spm' coordsystem

cfg = [];
template_seg = ft_volumesegment(cfg, template_mri); % segment the template mri

cfg = [];
cfg.method = 'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);   % create headmodel from segmented template
template_headmodel = ft_convert_units(template_headmodel, 'cm'); % convert to cm (CTF native coordinate system)

cfg = [];
cfg.grid.resolution = 1.5;
% cfg.grid.xgrid  = -40:1:40;
% cfg.grid.ygrid  = -40:1:40;
% cfg.grid.zgrid  = -40:1:40;
cfg.grid.tight = 'yes';
cfg.inwardshift = 0;                             % include region just outside head as 'inside', to avoid rim of nonactivation
cfg.headmodel = template_headmodel;
template_grid = ft_prepare_sourcemodel(cfg);        % create the template grid (coords of interest)

save('temp_grid','template_grid');
load temp_grid
%% 
% create headmodel from individual MRI

mripath = 'A:\controls\CTL01\nMR\T1';
individual_mri = ft_read_mri([mripath,'\T1.nii']);

% T1 = '.\nMR\T1\T1.nii';
% individual_mri = ft_read_mri(T1, 'dataformat', 'nifti_spm');

cfg = [];
cfg.method = 'flip';
individual_mri = ft_volumereslice(cfg, individual_mri);

% set fiducials

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
individual_mri = ft_volumerealign(cfg, individual_mri);

cfg = [];
cfg.output = 'brain';
individual_seg = ft_volumesegment(cfg, individual_mri);

 
% create the individual grid from mni grid
cfg = [];
cfg.grid.warpmni = 'yes';
cfg.grid.template = template_grid;
cfg.grid.nonlinear = 'yes';
cfg.mri = individual_mri;
individual_grid = ft_prepare_sourcemodel(cfg);
%% build the individual headmodel (needed for source analyses)

% create indidivual volume conductor model, for plotting

cfg = [];
cfg.method = 'singleshell';
individual_headmodel = ft_prepare_headmodel(cfg, individual_seg);
individual_headmodel = ft_convert_units(individual_headmodel, 'cm');

% sanity check - compare warped grid to headmodel built on individual anatomy

figure;
% subplot(1, 2, 1); 
ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none'); 
alpha 0.5; 
camlight; 
hold on; 
ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));

save('ana_sub1','individual_grid','individual_mri','individual_headmodel');

%%
dataset = 'A:\controls\CTL01\MEG\CTL01_verbs.ds';

figure,
ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.8); 
alpha 0.5; camlight; hold on; 
ft_plot_mesh(individual_grid.pos(individual_grid.inside, :)); 
sens = ft_read_sens(dataset); 
ft_plot_sens(sens);

%% Before
figure; hold on;
ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_mesh(individual_grid.pos(individual_grid.inside,:));

%%
atlas = ft_read_atlas('ROI_MNI_V4.nii'); 
atlas = ft_convert_units(atlas,'cm');
 
cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg, template_grid);

%% After
individual_grid2                 = individual_grid;
individual_grid2.inside          = false(individual_grid2.dim);
individual_grid2.inside(mask==1) = true;

figure; hold on;
ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_mesh(individual_grid2.pos(individual_grid2.inside,:));
 