clear; clc, close('all'); warning off


%% Initial settings
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
% indir = '/data/MEG/Clinical/MEG';
% %- Output dir 
% outdir = '/data/MEG/Clinical';

%-- Input dir
indir = '/data/MEG/Vahab/test_data';

%-- Output dir
outdir = '/data/MEG/Vahab/test_data/processed';


%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

switch task
    case 1
        % - Auditory definition naming
        tag = 'dfn';
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
        disp('5: Brainstorm source analysis using ft preprocessed');
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

%% analysis flag
flag.freq = 0;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis

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

%% All data
clear datafolder datafile
datafile1 = [];
d = rdir([indir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]);
d = rdir([indir,['/**/','/*',tag,'*/*raw_tsss.fif']]);

for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
end
datafile1 = vertcat(datafile1,datafile);
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
for i = 1:size(datafile1,1)
    
    datafile = datafile1{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(datafile)
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
    if flag.freq == 1
        Run_freq
    end
    
    %%
    cln_data_BAK = cln_data;
    
    %%
    cln_data = cln_data_BAK;
%     cfg                = [];
%     %     cfg.hpfilter       = 'yes';        % enable high-pass filtering
%     cfg.lpfilter       = 'yes';        % enable low-pass filtering
%     %     cfg.hpfreq       = 1;           % set up the frequency for high-pass filter
%     cfg.lpfreq         = 8;          % set up the frequency for low-pass filter
%     cln_data          = ft_preprocessing(cfg,cln_data);
    
    %     cfg             = [];
    % %     cfg.hpfilter    = 'yes';        % enable high-pass filtering
    %     cfg.lpfilter    = 'yes';        % enable low-pass filtering
    % %     cfg.hpfreq      = 12;           % set up the frequency for high-pass filter
    %     cfg.lpfreq      = 23;          % set up the frequency for low-pass filter
    %     cln_data        = ft_preprocessing(cfg,cln_data);
    
    %     cfg         = [];
    %     cfg.bsfilter = 'yes';
    %     cfg.bsfreq = [9 11]; % or whatever you deem appropriate
    %     cln_data   = ft_preprocessing(cfg,cln_data);
    
    %% Timelock
    if flag.time == 1
        switch task
            case 1
                %                 toi = [-0.5,0;0.3,0.8]; % Best of DN
                toi = [-0.3,0;1.1,1.8]; % Best of DN
            case 2
                toi = [0,0.3;0.6,1.2]; % Best for PN, left IFG
        end
        Run_time
    end
    
    %%
    %     cfg = [];
    %     cfg.savefile = [];
    %     cfg.saveflag = 2;
    %     cfg.lay  = lay;
    %     tfr = vy_tfr2(cfg, t_data);
    
    %% Grand Mean
    if flag.gave == 1
        Run_grandmean
    end
    %% Source analysis
    outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    %%
    switch analysis
        case 1
            %% Surface-based analysis
            vy_surfacebasedsource2
            vy_surfacebasedsource_dics
            
        case 2
            %% Volumetric-based analysis
            anatomy_check_flag = 2;
            Run_volumetric
        case 3
            %% SPM surface-based
            Run_spm   
        case 4
            %% BrainStorm export preprocessed ft-IC
            Run_bs
    end
    
    %%
    pause
    close all
    disp([datafile,' ,was completed'])
end





