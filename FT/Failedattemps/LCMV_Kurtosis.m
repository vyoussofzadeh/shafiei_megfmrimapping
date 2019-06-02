%%

cfg                  = [];
cfg.covariance       = 'yes';
test             = ft_timelockanalysis(cfg, ep_data.pst);

[u,s,v] = svd(test.cov);
figure;semilogy(diag(s),'o-');

%%
cfg                  = [];
cfg.method           = 'lcmv';
cfg.grid             = individual_grid;
cfg.headmodel        = individual_headmodel;
cfg.lcmv.keepfilter  = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.reducerank  = 2;
cfg.lcmv.lambda      = '5%';
cfg.lcmv.kappa       = 65;
cfg.lcmv.projectmom = 'yes';  %project dipole timeseries for each dipole in direction of maximal power (see below)
cfg.lcmv.kurtosis = 'yes';
source = ft_sourceanalysis(cfg, test);


%%
source.kurtosis = source.avg.kurtosis(source.inside); %get rid of NaNs which fall outside head
source.kurtosisdimord = 'pos';
source.pos     = template_grid.pos;
source.inside  = template_grid.inside;

cfg = [];
cfg.parameter = 'kurtosis';
source_interp = ft_sourceinterpolate(cfg, source, template_mri);


cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'kurtosis';
cfg.location     = param.loc;
ft_sourceplot(cfg,source_interp);

vy_mapvisualisation(source_interp,'kurtosis',0.5, []);


cfg = [];
cfg.funparameter = 'kurtosis';
cfg.method = 'ortho'; % orthogonal slices with crosshairs at peak (default anyway if not specified)
ft_sourceplot(cfg, source_interp);

cfg = [];
cfg.funparameter = 'kurtosis';
cfg.method = 'slice';  % plot slices
ft_sourceplot(cfg, source_interp);

