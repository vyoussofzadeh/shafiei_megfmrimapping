%% get data from SPM
spm12 = 'E:\My Matlab\SPM\spm12_3\spm12';
addpath(genpath(spm12));
% spm eeg

D = spm_eeg_load([path,megdata]);
% convert sensors and volume conduction model from SPM
volsens = spm_eeg_inv_get_vol_sens(D, 1, 'Head', 'inv', 'MEG');
vol1    = volsens.MEG.vol;
sens1   = volsens.MEG.sens;

disp(D.condlist)

% inv = D.inv;
% % mripath = inv{1, 1}.mesh.sMRI(1:end-2);
% mripath = inv{1, 1}.mesh.sMRI;
% fid = inv{1, 1}.mesh.fid.fid.pnt;
% mri_fids = fid(1:3,:);
 
% convert data from SPM
% raw      = D.ftraw(D.indchantype('MEGMAG'), D.indsample(-0.1):D.indsample(0.3), D.indtrial(D.condlist{1}, 'GOOD'));
% timelock = D.fttimelock(D.indchantype('MEGMAG'), D.indsample(-0.1):D.indsample(0.3), D.indtrial(D.condlist{1}, 'GOOD'));

% alternative method
raw      = spm2fieldtrip(D);

restoredefaultpath
%%
fieldtrip = 'E:\My Matlab\My codes\My GitHub\fieldtrip';
addpath(genpath(fieldtrip));

cfg = [];

cfg.chantype    = 'MEGMAG';
cfg.dftfilter = 'yes';
% cfg.hpfiltord = 3;
cfg.dftfreq = [60 120 180];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 1;
cfg.lpfreq = 40;
% cfg.channel = {'MEG'};
cfg.demean = 'yes';
cfg.baselinewindow = [-inf 0.0];
data_fix = ft_preprocessing(cfg,raw);

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
cfg.chantype    = 'MEGMAG';
comp = ft_componentanalysis(cfg, raw);

cfg = [];
c = 'no';
cfg.component = comp.label(51:end);
data_fix = ft_rejectcomponent(cfg, comp);


% cfg = [];
% cfg.method = 'pca';
% cfg.updatesens = 'no';
% cfg.chantype    = 'MEGMAG';
% comp = ft_componentanalysis(cfg, raw);
% 
% cfg = [];
% cfg.updatesens = 'no';
% cfg.component = comp.label(51:end);
% data_fix = ft_rejectcomponent(cfg, comp);

cfg = [];
cfg.chantype    = 'MEGMAG';
cfg.preproc.baselinewindow = [-inf 0];
% cfg.keeptrials = 'yes';
% cfg.preproc.bpfilter = 'yes';
% cfg.preproc.bpfreq = [1 30];
% cfg.preproc.hpfilter='yes';
% cfg.preproc.hpfreq=1;
% cfg.preproc.lpfilter='yes';
% cfg.preproc.lpfreq=40;
% cfg.preproc.demean = 'yes';
cfg.trials = find(data_fix.trialinfo==1); % Famous
timelock_fam = ft_timelockanalysis(cfg, data_fix);

cfg.trials = find(data_fix.trialinfo==2); % Unfamiliar
timelock_unf = ft_timelockanalysis(cfg, data_fix);

cfg.trials = find(data_fix.trialinfo==3); % Scrambled
timelock_scr = ft_timelockanalysis(cfg, data_fix);

cfg.trials = find(data_fix.trialinfo==4); % Face
timelock_fac = ft_timelockanalysis(cfg, data_fix);

cfg.trials = find(data_fix.trialinfo==5); % Faces - Scrambled
timelock_FS = ft_timelockanalysis(cfg, data_fix);

figure
cfg = [];
% cfg.layout = 'mpi_customized_acticap64.mat';
cfg.layout = 'neuromag306mag.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, timelock_fam, timelock_scr)

% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% difference = ft_math(cfg, timelock_fam, timelock_scr);
%%
cfg = [];
cfg.resamplefs = 300;
timelock_fam   = ft_resampledata(cfg, timelock_fam);
timelock_unf   = ft_resampledata(cfg, timelock_unf);
timelock_scr   = ft_resampledata(cfg, timelock_scr);
timelock_fac   = ft_resampledata(cfg, timelock_fac);
timelock_FS   = ft_resampledata(cfg, timelock_FS);
%%
cfg = [];
cfg.toilim = [-inf 0.0];
TL_Baseline = ft_redefinetrial(cfg, timelock_fam);
cfg.toilim = [0 0.5];
TL_post = ft_redefinetrial(cfg, timelock_fam);

% figure
% cfg = [];
% % cfg.layout = 'mpi_customized_acticap64.mat';
% cfg.layout = 'neuromag306mag.lay';
% cfg.interactive = 'yes';
% cfg.showoutline = 'yes';
% ft_multiplotER(cfg, TL_post)
% figure,
% ft_multiplotER(cfg, TL_Baseline)
%%
cfg = [];
cfg.covariance = 'yes'; % compute the covariance of the averaged ERF
% cfg.keeptrials = 'yes';
timelock_fam2 = ft_timelockanalysis(cfg, timelock_fam);
timelock_unf2 = ft_timelockanalysis(cfg, timelock_unf);
timelock_scr2 = ft_timelockanalysis(cfg, timelock_scr);
timelock_fac2 = ft_timelockanalysis(cfg, timelock_fac);
timelock_FS2 = ft_timelockanalysis(cfg, timelock_FS);

figure
plot(timelock_fam2.time, timelock_fam2.avg(71:376,:))

figure,
cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, timelock_fam2)

% figure,
% cfg = [];
% cfg.layout = 'neuromag306mag.lay';
% cfg.interactive = 'yes';
% cfg.showoutline = 'yes';
% ft_multiplotER(cfg, timelock_FS2)

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'no';
datapow          = ft_freqanalysis(cfg, timelock_fam2);

cfg        = [];
cfg.layout = 'neuromag306mag.lay';
cfg.xlim   = [10 20];
figure,ft_topoplotER(cfg, datapow);
