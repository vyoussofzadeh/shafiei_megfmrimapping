% BS2FT
pialdir = fullfile(datadir,subj,'brainstorm_db/anat',subj,'tess_cortex_pial_low.mat');
sourcemodel = ft_read_headshape(pialdir);
%     figure
%     ft_plot_mesh(sourcemodel, 'vertexcolor', sourcemodel.curv)

[individual_headmodel, individual_grid] = vy_bs2ft_headmodel(t_data, sourcemodel);

figure;
ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
view ([-10 40 10])

%     hsfile = datafile; % headshape
%     headshape = ft_read_headshape(hsfile);
%     headshape = ft_convert_units(headshape,'mm');
%     ft_plot_headshape(headshape);

%
method = 'lcmv';
cfg                  = [];
cfg.method           = 'lcmv';
cfg.grid             = individual_grid; % leadfield, which has the grid information
cfg.headmodel        = individual_headmodel; % volume conduction model (headmodel)
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.lambda      = '10%';
cfg.lcmv.keepfilter  = 'yes';
sourceAll = ft_sourceanalysis(cfg, t_data.app);
cfg.grid.filter = sourceAll.avg.filter;
s_data.bsl = ft_sourceanalysis(cfg, t_data.bsl);
s_data.pst = ft_sourceanalysis(cfg, t_data.pst);

s_data2.bsl  = ft_sourcedescriptives([], s_data.bsl); % to get the neural-activity-index
s_data2.pst  = ft_sourcedescriptives([], s_data.pst); % to get the neural-activity-index

%
cfg = [];
cfg.parameter = 'pow';
%     cfg.operation = '(x1-x2)/(x1+x2)';
cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
source_ratio = ft_math(cfg,s_data.pst,s_data.bsl);
%     source_ratio.pow(source_ratio.pow>0)=0;
%     source_ratio.pow = abs(source_ratio.pow);

figure
m = source_ratio.pow;
bnd.pnt = individual_grid.pos;
bnd.tri = individual_grid.tri;
ft_plot_mesh(bnd, 'vertexcolor', m);
colorbar