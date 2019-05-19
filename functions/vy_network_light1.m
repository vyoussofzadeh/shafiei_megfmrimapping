function vy_network_light1(cfg_main, t_data)

% check if lcmv exists!
% if ~exist('s_data2','var')
[s_data, s_data2] = vy_source(t_data, cfg_main.grid, cfg_main.headmodel);
% end
%%
% using older ver of ft for network analysis
restoredefaultpath
addpath(genpath(cfg_main.p.ft_old));
addpath(genpath(cfg_main.p.hcp_path));
addpath(genpath([cfg_main.p.cd_org,'/functions']));
addpath(genpath([cfg_main.p.cd_org,'/Data_file']));

%%
mtd = 'plv';
mtd_par = 'plvspctrm';
gtm = 'eigenvector_cent';

conn_par = [];
conn_par.method   = mtd;
conn_par.idx      = mtd_par;
conn_par.complex  = [];

net_par.gtm       = gtm ; in2 = 2;
net_par.threshold = 0.7;
[source_conn_bsl, network_bsl] = vy_conn(s_data2.bsl,conn_par,net_par);
[source_conn_pst, network_pst] = vy_conn(s_data2.pst,conn_par,net_par);

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
network_diff_lcmv = ft_math(cfg,network_pst,network_bsl);
network_diff_lcmv.pos     = cfg_main.template_grid.pos;
network_diff_lcmv.dim     = cfg_main.template_grid.dim;
network_diff_lcmv.inside  = cfg_main.template_grid.inside;
network_diff_lcmv.(gtm)   = zscore(network_diff_lcmv.(gtm));
network_diff_lcmv.eigenvector_cent(network_diff_lcmv.eigenvector_cent<0)=0;

%% revert to the newer ft!
restoredefaultpath
addpath(genpath(cfg_main.p.ft_path));
addpath(genpath(cfg_main.p.hcp_path));
addpath(genpath([cfg_main.p.cd_org,'/functions']));
addpath(genpath([cfg_main.p.cd_org,'/Data_file']));

%%
savepath = fullfile(cfg_main.outputdir);
if exist(savepath, 'file') == 0, mkdir(savepath), end
savedata = fullfile(savepath,['n_',cfg_main.subj,'.mat']);
save(savedata, 'network_diff_lcmv', '-v7.3');
mtd = 'network_evc';


% oldmask = network_diff_lcmv.(gtm);
% network_diff_lcmv.mask = (oldmask - min(oldmask(:))) ./ (max(oldmask(:)) - min(oldmask(:))); %
% % set the new mask to range between 0 and 1
% network_diff_lcmv.mask(oldmask < 0.9) = 0; % mask out non-significant voxels

param = [];
param.mask = gtm;
param.loc = 'max';
network_int_lcmv = vy_source_plot(network_diff_lcmv,cfg_main.template_mri,param,2);
savefig = fullfile(savepath,[mtd,'_1_',cfg_main.subj]);
hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);

% clear savep
% savep{1} = fullfile(savepath,[mtd,'_2_',cfg_main.subj]);
% savep{2} = fullfile(savepath,[mtd,'_3_',cfg_main.subj]);
% vy_mapvisualisation(network_int_lcmv,gtm,0.6, savep);
vy_mapvisualisation(network_int_lcmv,gtm,0.6, []);


% savenii = fullfile(savepath,['n_',subj,'.nii']);
% vy_savenifti(network_int_lcmv, gtm, savenii);

%%
% cfg = [];
% cfg.funparameter = param.mask;
% cfg.method = 'ortho';  % plot slices
% ft_sourceplot(cfg, network_diff_lcmv);


%%
% network_int1 = network_int_lcmv;
% thre = 0.5;
% network_int1.eigenvector_cent(network_int1.eigenvector_cent < thre*max(network_int1.eigenvector_cent(:)))=0;
% % network_int1.eigenvector_cent(network_int1.eigenvector_cent > 0) = 1;
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = gtm;
% % cfg.filename  = './output';
% cfg.filename  = fullfile(outputdir_net,'surf.nii');
% ft_volumewrite(cfg, network_int1)
%
% restoredefaultpath
% addpath(genpath([cd_org,'/functions']));
% % close all
% addpath(connpath);
% addpath(genpath(spm_path))
% h = get(0, 'Children');
% if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
%     conn
% end
%
% filenameVOL = [];
% FSfolder = fullfile(connpath,'/utils/surf');
% sphplots = [];
% connplots = [];
% facealpha = 1;
% position = [-1 0  0];
% conn_mesh_display(fullfile(outputdir_net,'surf.nii'),[],FSfolder);

%% revert to the newer ft!
% restoredefaultpath
% addpath(genpath(cfg_main.p.ft_path));
% addpath(genpath(cfg_main.p.hcp_path));
% addpath(genpath([cfg_main.p.cd_org,'/functions']));
% addpath(genpath([cfg_main.p.cd_org,'/Data_file']));



%% parcellation - aal (132 rois)
% network_diff_lcmv.eigenvector_cent = (network_diff_lcmv.eigenvector_cent)./max(network_diff_lcmv.eigenvector_cent);
% [~, data_intpar, coor] = vy_parcellate(network_diff_lcmv, atlas, gtm);
% data_intpar.eigenvector_centdimord = 'chan';
%

%% Conn
% conn_par.conn_thre = 0.95;
% conn_ratio = vy_connvis(source_conn_pst,source_conn_bsl,conn_par, individual_headmodel, network_diff_lcmv);
% view([156,47]);

%
%% ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, gtm);
% disp(ROI_sel)
% savepath = fullfile(outputdir_net,['n_ROIs_',subj]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);


end