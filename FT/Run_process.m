clear; clc, close('all'); warning off

%%
restoredefaultpath

% Initial settings, tools and data directory
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));

cd_org = cd;
cd_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
cd (cd_tools)

vy_init

datadir = '/data/MEG/Clinical/MEG';
cd(cd_org)

addpath(genpath(cd_org));

outdir = '/data/MEG/Clinical';
% cd(outdir)

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

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
disp('1: Surface-based')
disp('2: Volumetric');
analysis = input('Eneter the analysis: ');

switch analysis
    case 1
        mtag = 'source_surf';
    case 2
        % end
        disp('1: LCMV source')
        disp('2: Network/Connectvity - Broadband');
        disp('3: DICS Source, Beta');
        disp('4: SPM source analysis (surface + BF)');
        method = input('Eneter the method: ');
        switch method
            case 1
                mtag = 'lcmv';
            case 2
                mtag = 'conn';
            case 3
                mtag = 'dics';
        end
        %-
        disp('1: Low-res grid')
        disp('2: High-res grid')
        meshgrid = input('Eneter the mesh grid: ');
end

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: older');
year = input('Year data was acquired: ');

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

%% all data
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);

%% Per year
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

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%%
for i = 2:size(datafile1,1)
    %%
    datafile = datafile1{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(['subj:',subj])
    disp(['Date:',Date])
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    
    %%
    if year==4
        yttag = 'older';
    else
        yttag = ytag{1};
    end
    outd.sub = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    
    %% Preprocesssing
    Run_preprocess
    
    %% FFT & TFR
    Run_freq
    
    %% Timelock
%     toi = [-0.3,0;0.7,1.8];
    toi = [-0.3,0;0.7,1.8];
    Run_time
    
    %% Grand Mean
    Run_grandmean
    
    %%
    outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
%%
    switch analysis
        case 1
            %% Surface-based analysis
            outd.sourcesuf = fullfile(outd.sub,mtag);
            cfg = [];
            cfg.task = tag;
            cfg.outputdir = outd.sourcesuf;
            cfg.subj = subj;
            cfg.data = t_data.pst;
            cfg.datadir = datadir;
            cfg.outputmridir = outputmridir;
            cfg.peaksel = 4;
            vy_surfacebasedsource(cfg)
            
        case 2
            %%
            Run_Volumetric
    end
    
    %%
    pause
    close all
    disp([datafile,' ,was completed'])
end





