function freq = vy_fft(cfg_mian, data)

%%
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'pow';
cfg.pad          = 'nextpow2';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_mian.foilim;
cfg.tapsmofrq    = 4;
cfg.taper        = 'hanning';
cfg.pad          = 4;
freq             = ft_freqanalysis(cfg, data);

if cfg_mian.plotflag ==1
    psd = squeeze(mean(mean(abs(freq.fourierspctrm),2),1));
    ff = linspace(1, cfg.foilim(2), length(psd));
    figure,plot(ff,psd)
    xlabel('Hz'); ylabel('psd')
    max(psd)
end


if cfg_mian.saveflag ==1
    hcp_write_figure([cfg_mian.savefile,'.png'], gcf, 'resolution', 300); 
end