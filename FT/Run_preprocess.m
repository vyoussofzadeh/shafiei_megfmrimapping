
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2   
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
else
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outd.sub,['f_',subj,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        event = ft_read_event(datafile);
        clear val
        for i=1:length(event)
            val(i) = event(i).value;
        end
        val1 = unique(val);
        %
        disp('preprocessing ...');
        cfg = [];
        cfg.eventid = min(val1);
        cfg.epochtype = event(1).type;
        cfg.datafile  = datafile;
        cfg.hpfreq = 0.1;
        cfg.lpfreq = 40;
        [f_data, ecg_data] = vy_preprocess(cfg);
        disp('filtering was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
    %% Bad trials & channels (automated)
    savepath = fullfile(outd.sub,['a_',subj,'.mat']);
    cfg = [];
    cfg.pflag = 1; % yes:1, No:2
    cfg.saveflag = 1; % yes:1, No:2
    cfg.savepath = savepath;
    cfg.latency = [-200,900];
    [r_data,report] = vy_artifactreject(cfg, f_data);
    % disp('Bad data rejection was completed');
    
    %% Bad trials & channels (Manuual)
    % clear r_data
    % cfg = [];
    % cfg.metric = 'zvalue';  % use by default zvalue method
    % cfg.latency = [-200,900];
    % cfg.layout   = lay;   % this allows for plotting individual trials
    % r_data   = ft_rejectvisual(cfg, f_data);
    
    %% Inspecting bad data
    % cfg = [];
    % cfg.viewmode = 'vertical';
    % cfg.continuous = 'no';
    % cfg.trials     = report.btrl;
    % cfg.channel   = report.bchan;
    % ft_databrowser(cfg,f_data);
    
    %% ICA cleaning
    cfg = [];
    cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.lay = lay;
    cfg.n   = 20;
    cln_data = vy_ica_cleaning_light(cfg, r_data);
    disp('ICA cleaning was completed');
    
end