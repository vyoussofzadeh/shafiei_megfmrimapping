%%
% vol1 = ft_convert_units(vol1, 'mm');
% sens1 = ft_convert_units(sens1, 'mm');
load(fid);

%%
% dataset_raw = output.dataID.dataset_raw;
% dataset = output.dataID.dataset;
% data_meg = output.preprocess.data;

%% start from scratch data in FieldTrip

% subj = 15;
% prefix = sprintf('Sub%02d', subj);
% load([prefix '_raw']);  % this is called "data" rather than "raw"

% sens = data.grad;

% load the original MRI
mri_orig = ft_read_mri(mripath);

% load the positions of the anatomical fiducials (as provided by Rik)
% load('/home/common/matlab/fieldtrip/data/test/original/rikhenson/Sub15/T1/mri_fids.mat');

% headshape = ft_read_headshape(dataset_raw);
% headshape = ft_convert_units(headshape, 'mm');

% the MRI is neither expressed in MNI, nor in Neuromag coordinates
ft_determine_coordsys(mri_orig, 'interactive', 'no');
hold on; % add the subsequent objects to the same figure
% ft_plot_headshape(headshape);
plot3(mri_fids(1,1), mri_fids(1,2), mri_fids(1,3), 'm*');
plot3(mri_fids(2,1), mri_fids(2,2), mri_fids(2,3), 'm*');
plot3(mri_fids(3,1), mri_fids(3,2), mri_fids(3,3), 'm*');
rotate3d

%% validate the positions of the fiducials that were provided by Rik

% cfg = [];
% cfg.location = mri_fids(1,:);
% ft_sourceplot(cfg, mri_orig);
% 
% cfg = [];
% cfg.location = mri_fids(2,:);
% ft_sourceplot(cfg, mri_orig);
% 
% cfg = [];
% cfg.location = mri_fids(3,:);
% ft_sourceplot(cfg, mri_orig);

%%

% the location of fiducials is expressed in original MRI coordinates
% ft_volumerealign needs them in voxel coordinates
vox_fids = ft_warp_apply(inv(mri_orig.transform), mri_fids);

cfg = [];
cfg.fiducial.nas = vox_fids(1,:);
cfg.fiducial.lpa = vox_fids(2,:);
cfg.fiducial.rpa = vox_fids(3,:);
cfg.coordsys = 'neuromag';
mri_realigned = ft_volumerealign(cfg, mri_orig);

% save mri_realigned mri_realigned

% check that the MRI is consistent after realignment
% ft_determine_coordsys(mri_realigned, 'interactive', 'no');
% hold on; % add the subsequent objects to the figure
% ft_plot_headshape(headshape);

%%
% cfg = [];
% cfg.output = {'brain' 'scalp' 'skull'};
% seg = ft_volumesegment(cfg, mri_realigned);
% 
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.numvertices = 2000;
% cfg.tissue = 'brain';
% brain = ft_prepare_mesh(cfg, seg);
% cfg.tissue = 'skull';
% skull = ft_prepare_mesh(cfg, seg);
% cfg.tissue = 'scalp';
% scalp = ft_prepare_mesh(cfg, seg);

%%
% cfg = [];
% cfg.method = 'singleshell';
% vol = ft_prepare_headmodel(cfg, brain);

% save vol vol
% save sens sens

% ft_determine_coordsys(mri_realigned, 'interactive', 'no')
% hold on; % add the subsequent objects to the same figure
% ft_plot_headshape(headshape);
% ft_plot_vol(ft_convert_units(vol, 'mm'));

% figure
% hold on; % add the subsequent objects to the same figure
% % ft_plot_headshape(headshape);
% ft_plot_sens(sens1);
% 
% % ft_plot_sens(ft_convert_units(sens1, 'mm'), 'coil', 'yes', 'coildiameter', 10);
% ft_plot_vol(vol1);
% 
% figure
% ft_plot_vol(ft_convert_units(vol1,  'mm'), 'facecolor', 'r'); % FT
% alpha 0.5
