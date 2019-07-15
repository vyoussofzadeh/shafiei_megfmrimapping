f = 23;
tapsmofrq = 4;
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;
if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = sens;
%%

cfg = [];
cfg.method = 'dics';
cfg.dics.lambda = '0%';
cfg.frequency    = f;
cfg.grid = individual_grid;
cfg.headmodel = individual_headmodel;
% cfg.dics.keepfilter  = 'yes';
% cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.dics.projectnoise = 'yes';
sourceavg = ft_sourceanalysis(cfg, f_data.pst);

%%
sourceNAI = sourceavg;
sourceNAI.avg.pow = sourceavg.avg.pow ./ sourceavg.avg.noise;
sourceNAI.pos     = template_grid.pos;
sourceNAI.dim     = template_grid.dim;
sourceNAI.inside  = template_grid.inside;

%%

cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
source_int_dics = vy_source_plot(cfg, sourceNAI);