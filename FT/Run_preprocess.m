

%% Filteting, Event reading, Artifact rejecrtion
savepath = fullfile(outputdir,['f_',subj,'.mat']);
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
    [f_data, ecg_data] = vy_preprocess(cfg);
    disp('preprocessing was completed');
    save(savepath, 'f_data', '-v7.3');
end

%% Visual artifacts
savepath = fullfile(outputdir,['a_',subj,'.mat']);
cfg = [];
cfg.pflag = 1; % yes:1, No:2
cfg.saveflag = 1; % yes:1, No:2
cfg.savepath = savepath;
r_data = vy_artifactreject(cfg, f_data);

%% ICA cleaning
cfg = [];
cfg.savepath = fullfile(outputdir,['ic_',subj,'.mat']);
cfg.saveflag = 1;
cfg.lay = lay;
cfg.n   = 20;
ic_data = vy_ica_cleaning_light(cfg, r_data);







