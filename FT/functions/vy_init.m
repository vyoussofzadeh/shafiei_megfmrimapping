
%- fieldtrip
ft_path = fullfile(cd_tools,'/ft packages/fieldtrip-master');
ft_path_2019 = fullfile(cd_tools,'/ft packages/fieldtrip_20190419');

addpath(ft_path);
ft_defaults

ft_old = fullfile(cd_tools,'ft packages/ft_old');
%- colormap
addpath(fullfile(ft_path, 'external/brewermap'));
%-fastICA
addpath(ft_path,'/external/fastica');

%- HCP
hcp_path = fullfile(cd_tools,'/megconnectome-3.0');
addpath(genpath(hcp_path));

%- atlas
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));

%-Grid template
load temp_grid_8mm % from, vy_warping()
% template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %
template_mri = ft_read_mri(fullfile(ft_path,'template/anatomy','single_subj_T1.nii')); %

%-CONN
connpath = fullfile(cd_tools,'/tools/Conn/conn');
% addpath(connpath);

%-SPM
% spm_path = '/opt/matlab_toolboxes/spm12';
spm_path = fullfile(cd_tools,'SPM/spm12');
% addpath(spm_path);

%-SPM-beamformer
spmbf_path = fullfile(cd_tools,'Beamforming');

%%
% % ft_path = 'F:/My Matlab/My codes/My GitHub/fieldtrip_041718/fieldtrip-master';
% ft_path = 'F:/My Matlab/My codes/My GitHub/fieldtrip_packs/fieldtrip-20180809';
% 
% addpath(genpath(fullfile(cd_org,'functions')));
% addpath(genpath(fullfile(cd_org,'Data_file')));
% 
% 
% % ft_path = 'F:/My Matlab/My codes/My GitHub/fieldtrip_new/fieldtrip-20170321';
% addpath(ft_path);
% ft_defaults % this loads the rest of the defaults;
% 
% hcp_path = 'F:/My Matlab/MEG/HCP/megconnectome-3.0';
% addpath(genpath(hcp_path));
% 
% atlas_path = fullfile(ft_path,'template','atlas');
% atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));
% 
% headmodel_file = 'standard_bem.mat';
% % load temp_grid_8mm % from, vy_warping()
% template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %
% 
% spm_path = 'F:/My Matlab/SPM/spm12_4/spm12/spm12';
% addpath(genpath(spm_path))
% 
% %- colormap
% addpath(fullfile(ft_path, 'external/brewermap'));
% 
% %-ica
% % addpath(ft_path,'external/eeglab');
% 
% %-fastICA
% addpath(ft_path,'/external/fastica');
% 
% % - trl functions
% % addpath(genpath('F:/My Matlab/Persons/ABF/MEG_processing_scripts/TrialFun'));
