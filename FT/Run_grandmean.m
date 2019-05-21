a_data = vy_ave(cln_data);
savepath = fullfile(outd.sub,'Timelock');
if exist(savepath, 'file') == 0, mkdir(savepath), end
cfg = [];
cfg.savefile = fullfile(savepath,['gmean_',subj,'.mat']);
cfg.saveflag = 1;
cfg.lay  = lay;
vy_ave_plot(cfg, a_data);