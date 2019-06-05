function source_int = vy_source_plot(data,template_mri,param,n)

cfg = [];
cfg.parameter = param.mask;
% cfg.interpmethod = 'sphere_avg';
cfg.interpmethod = 'smudge';
cfg.coordsys     = 'mni';
source_int = ft_sourceinterpolate(cfg, data, template_mri);

if n==1
    cfg = [];
    source_int = ft_volumenormalise(cfg, source_int);
end

cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = param.mask;
cfg.location     = param.loc;
cfg.funcolormap =  brewermap(256, '*RdYlBu');
ft_sourceplot(cfg,source_int);

%%
source_int1 = source_int;
tmp = abs(source_int1.(param.mask));
tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
source_int1.(param.mask) = tmp;

% plot multiple 2D axial slices
cfg = [];
% cfg.method        = 'ortho';
cfg.method        = 'slice';
cfg.funparameter  = param.mask;
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.1 1];
cfg.opacitylim    = [0.1 1];
cfg.nslices       = 12;
cfg.slicerange    = [20,50];
cfg.slicedim      = [3];
cfg.opacitymap    = 'rampup';
cfg.funcolormap =  brewermap(256, '*RdYlBu');
ft_sourceplot(cfg, source_int1);


%%
% A_range = [0 5];
% % Use the T-statistics to create an alpha map (which must be in [0,1])
% alphamap = abs(source_int.(param.mask));
% alphamap(alphamap > A_range(2)) = A_range(2);
% alphamap(alphamap < A_range(1)) = 0;
% alphamap = alphamap/A_range(2);
% 
% 
% 
% cfg              = [];
% % cfg.method       = 'ortho';
% cfg.method        = 'slice';
% cfg.funparameter = param.mask;
% cfg.funcolorlim   = [2 3];
% cfg.opacitylim    = [2 3];
% cfg.opacitymap    = 'rampup';
% cfg.location     = param.loc;
% ft_sourceplot(cfg,source_int);
% 
% H.AlphaData = alphamap;
% 
% 
% H = get(gca,'Children');
% % Adjust the alpha values of the overlay 
% set(H, 'AlphaData', alphamap);



% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funcolormap   = 'jet';
% cfg.location      = 'max';
% cfg.funparameter  = 'pow';
% cfg.funcolormap = flipud(brewermap(64,'RdBu'));
% cfg.crosshair = 'no';
% ft_sourceplot(cfg, source_int);
% set(gca,'color','none')
