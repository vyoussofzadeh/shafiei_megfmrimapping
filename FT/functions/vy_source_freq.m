function [s_data, s_data2] = vy_source_freq(data, grid, vol, mtd)


switch mtd
    
    case 'dics'
        cfg              = [];
        cfg.method       = 'dics';
        cfg.frequency    = data.app.freq;
        cfg.grid         = grid;
        cfg.headmodel    = vol;
        cfg.dics.lambda       = '5%';
        cfg.dics.keepfilter   = 'yes';
        %         cfg.dics.realfilter   = 'yes';
        cfg.dics.fixedori     = 'yes';
        sourceAll = ft_sourceanalysis(cfg, data.app);
        cfg.grid.filter = sourceAll.avg.filter;
        s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst = ft_sourceanalysis(cfg, data.pst);
        
    case 'dics_stat'
        cfg              = [];
        cfg.method       = 'dics';
        cfg.frequency    = data.app.freq;
        %         cfg.grid         = grid;
        cfg.grid.pos = grid.pos(grid.inside,:);
        cfg.headmodel    = vol;
        cfg.dics.lambda       = '5%';
        cfg.dics.keepfilter   = 'yes';
        cfg.dics.fixedori     = 'yes';
        sourceAll = ft_sourceanalysis(cfg, data.app);
%         cfg.grid.filter = sourceAll.avg.filter;
%         cfg.rawtrial = 'yes';
%         cfg.vol = vol;
%         s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
%         s_data.pst = ft_sourceanalysis(cfg, data.pst);
                
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '5%';
        cfg.frequency    = data.app.freq;
        % cfg.grid = grid;
        cfg.grid.pos = grid.pos(grid.inside,:);
        cfg.headmodel = vol;
        cfg.dics.keepfilter = 'yes';
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.channel = data.bsl.label;
        sourceavg = ft_sourceanalysis(cfg, data.app);
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '5%';
        cfg.grid.pos = grid.pos(grid.inside,:);
        cfg.grid.filter = sourceavg.avg.filter;
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.rawtrial = 'yes';
        cfg.headmodel = vol;
        s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, data.pst);
        
    case 'pcc'
        
        cfg              = [];
        cfg.method       = 'pcc';
        cfg.frequency    = data.app.freq;
        cfg.grid         = grid;
        cfg.headmodel    = vol;
        cfg.pcc.lambda       = '5%';
        cfg.pcc.keepfilter   = 'yes';
        cfg.pcc.projectnoise  = 'yes';
        cfg.pcc.fixedori     = 'yes';
        sourceAll = ft_sourceanalysis(cfg, data.app);
        cfg.grid.filter = sourceAll.avg.filter;
        s_data.bsl = ft_sourceanalysis(cfg, data.bsl);
        s_data.pst = ft_sourceanalysis(cfg, data.pst);
end

cfg = [];
% cfg.projectmom  = 'yes';
s_data2.bsl  = ft_sourcedescriptives(cfg, s_data.bsl); % to get the neural-activity-index
s_data2.pst  = ft_sourcedescriptives(cfg, s_data.pst); % to get the neural-activity-index




