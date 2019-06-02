function comp = vy_ica(data,lay, n)

cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = n;       % specify the component(s) that should be plotted
comp           = ft_componentanalysis(cfg, data);


%%
cfg           = [];
cfg.component = 1:n;       % specify the component(s) that should be plotted
cfg.layout    = lay;
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)
colormap(brewermap(256, '*RdYlBu'));


cfg = [];
cfg.viewmode = 'component';
cfg.layout = lay;
ft_databrowser(cfg, comp);
colormap(brewermap(256, '*RdYlBu'));
% set(gcf, 'Position', [600   600   700   500]);

%%
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';%compute the power spectrum in all ICs
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;
freq = ft_freqanalysis(cfg, comp);


%%
nby1 = 5; nby2 = 4;

Nfigs = ceil(size(comp.topo,1)/n);
tot = Nfigs*n;

rptvect = 1:size(comp.topo,1);
rptvect = padarray(rptvect, [0 tot-size(comp.topo,1)], 0,'post');
rptvect = reshape(rptvect,n,Nfigs)';

figure
for r=1:n
    cfg=[];
    cfg.channel = rptvect(:,r);
    subplot(nby1,nby2,r);set(gca,'color','none');
    ft_singleplotER(cfg,freq);
end
colormap(brewermap(256, '*RdYlBu'));
set(gcf, 'Position', [800   600   800   500]);

end
