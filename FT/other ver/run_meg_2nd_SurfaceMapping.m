clear; clc, close('all')
%%
addpath(genpath('.\functions'));
vy_init
addpath('.\Group results');
addpath(connpath);
addpath(spm_path);

%%

%-fmri
vol_fmri_CRM = ft_read_mri('fmri_CRM.nii');
Opt = [];
Opt.savenii = 0;
Opt.savefig = 0;
Opt.savename = '.\Group results\fmri_CRM';
vy_surfce_vis2(vol_fmri_CRM,'.\Group results\fmri_CRM.nii', Opt);

vol_fmri_VGA = ft_read_mri('fmri_VGA.nii');
Opt = [];
Opt.savenii = 0;
Opt.savefig = 0;
Opt.savename = '.\Group results\fmri_VGA';
vy_surfce_vis2(vol_fmri_VGA,'.\Group results\fmri_VGA.nii', Opt);

%-meg
load ('meg_CRM.mat');

s_vol_meg = vy_vol_thresh(s_vol_meg,projthresh_meg); % abs
s_vol_meg_group = vy_vol_thresh_posneg(s_vol_meg,.001, 1); % keep positive effects only!

Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = '.\Group results\meg_CRM';
vy_surfce_vis2(s_vol_meg_group,Opt.savename, Opt);


projthresh = 0.6; 
vol_meg_CRM = ft_read_mri('meg_CRM.nii');
vol_meg_CRM = vy_vol_thresh(vol_meg_CRM,projthresh); % abs
Opt = [];
Opt.savenii = 1;
Opt.savefig = 1;
Opt.savename = '.\Group results\meg_CRM';
vy_surfce_vis2(vol_meg_CRM,'.\Group results\meg_CRM.nii', Opt);

vol_meg_VGA = ft_read_mri('meg_VGA.nii');
Opt = [];
Opt.savenii = 0;
Opt.savefig = 0;
Opt.savename = '.\Group results\meg_VGA';
vy_surfce_vis2(vol_meg_VGA,'.\Group results\meg_VGA.nii', Opt);


%% Masking combined meg-fmri
%-inspired by, https://www.nitrc.org/forum/forum.php?thread_id=7218&forum_id=1144
savemask = fullfile(savedir,['meg_net_fmri_',tsk,'.nii']);
spm_imcalc({vol_name_meg, vol_name_fmri},savemask, 'i1.*(i1>0)-i2.*(i2>0)');

Opt = [];
Opt.savenii = 2;
Opt.savefig = 1;
Opt.savename = fullfile(savedir,['meg_net_fmri_',tsk]);
vy_surfce_vis2([],savemask, Opt);
