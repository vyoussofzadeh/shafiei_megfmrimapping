function vy_mapvisualisation(input,funparameter,thre,savepath)

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = funparameter;
cfg.funcolorlim = 'maxabs';
cfg.opacitymap = 'rampup';
cfg.crosshair = 'no';
cfg.camlight       = 'no';
cfg.funcolormap =  brewermap(256, '*RdYlBu');
cfg.projthresh     = thre;

cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.surfinflated   = 'surface_inflated_right.mat';
ft_sourceplot(cfg, input);
% view ([-70 20 50])
% light ('Position',[-70 20 50]);
view([90 0]); camlight; material dull; 
% title('test');

if ~isempty(savepath)==1
    hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);
end

ft_sourceplot(cfg, input);
% view ([70 20 50])
% light ('Position',[70 20 50])
view([-90 0]); camlight; material dull; 

if ~isempty(savepath)==1
    hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);
end



