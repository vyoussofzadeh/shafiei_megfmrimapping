clear; clc, close('all'); warning off

%%
restoredefaultpath

%% Initial settings, tools and data directory
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));

cd_org = cd;
cd_tools = '/data/MEG/Vahab/Scripts/Vahab/Scripts/tools';
cd (cd_tools)

vy_init

datadir = '/data/MEG/Clinical/MEG';
cd(cd_org)

outdir = '/data/MEG/Clinical';
% cd(outdir)

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

%%
switch task
    case 1
        % - Auditory definition naming
        tag = 'DFN';
        Evnt_IDs = 1; % questions
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
end

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: older');

year = input('Year data acquired: ');

clear ytag;
switch year
    case 1
        ytag = {'19'};
    case 2
        ytag = {'18'};
    case 3
        ytag = {'17'};
    case 4
        ytag = {'up','11','12','13','14','15','16'};
end

%%
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
clear datafolder datafile
datafile1 = [];
for j=1:numel(ytag)
    ytag1 = ytag{1,j};
    d = rdir([datadir,['/**/',ytag1,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
    end
    datafile1 = vertcat(datafile1,datafile);
end
datafile1 = datafile1';
disp(datafile1)

%%
epoch_type = 'STI101';
% run = '1';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%%
for i = 1:size(datafile1,1)
    
    datafile = datafile1{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(['subj:',subj])
    disp(['Date:',Date])
    
    %%
    if year==4
        yttag = 'older';
    else
        yttag = ytag{1};
    end
    
    outputdir = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outputdir, 'file') == 0
        mkdir(outputdir);   %create the directory
    end
    cd(outputdir)
    disp(['outputdir:',outputdir])
    
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
    savepath = fullfile(outputdir,['r_',subj,'.mat']);
    cfg = [];
    cfg.pflag = 1;
    cfg.saveflag = 2;
    cfg.savepath = savepath;
    r_data = vy_artifactreject(cfg, f_data);
    
    %% ICA cleaning
    savepath = fullfile(outputdir,['ica_',subj,'.mat']);
    saveflag = 1;
    ic_data = vy_ica_cleaning(r_data, lay, savepath, saveflag);
        
end


%%
%     cfg = [];
%     cfg.artfctdef.zvalue.cutoff = 7;
%     cfg.artfctdef.zvalue.channel = f_data.label;
%     cfg.continuous  = 'no';
%     Z_artifact = ft_artifact_zvalue(cfg, f_data);
%
%     cfg = [];
%     cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
%     cfg.artfctdef.zvalue.artifact = Z_artifact;
%     f_data1 = ft_rejectartifact(cfg,f_data);

%
%     cfg = [];
%     cfg.trials = 'all';
%     cfg.metric = 'kurtosis';
%     f_data2 = ft_rejectvisual(cfg, f_data1);






