clear; clc, close('all'); warning off

%%
restoredefaultpath

% Initial settings, tools and data directory
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
disp('1: Surface-based')
disp('2: Volumetric');
analysis = input('Eneter the analysis: ');

%%
disp('1: Low-res grid')
disp('2: High-res grid')
meshgrid = input('Eneter the mesh grid: ');

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
        cfg.hpfreq = 1;
        cfg.lpfreq = 28;
        [f_data, ecg_data] = vy_preprocess(cfg);
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
    %% Bad trials & channels
    savepath = fullfile(outputdir,['a_',subj,'.mat']);
    cfg = [];
    cfg.pflag = 1; % yes:1, No:2
    cfg.saveflag = 1; % yes:1, No:2
    cfg.savepath = savepath;
    cfg.latency = [-400,1200];
    [r_data,report] = vy_artifactreject(cfg, f_data);
    
    %% Inspecting bad data
%     cfg = [];
%     cfg.viewmode = 'vertical';
%     cfg.continuous = 'no';
%     cfg.trials     = report.btrl;
%     cfg.channel   = report.bchan;
%     ft_databrowser(cfg,f_data);
    
    %% ICA cleaning
    cfg = [];
    cfg.savepath = fullfile(outputdir,['ic_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.lay = lay;
    cfg.n   = 20;
    cln_data = vy_ica_cleaning_light(cfg, r_data);
    
    %% freq analysis (fft)
    savepath = fullfile(outputdir,'Freq');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    cfg = [];
    cfg.savefile = fullfile(savepath,['fft_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.foilim = [2 40];
    cfg.plotflag  = 1;
    vy_fft(cfg, cln_data);
    
    %% Time-freq analysis (tfr)
    savepath = fullfile(outputdir,'Freq');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    cfg = [];
    cfg.savefile = fullfile(savepath,['tfr_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.lay  = lay;
    tfr = vy_tfr(cfg, cln_data);
    
    % tfr plotting
    cfg = [];
    cfg.savepath = fullfile(outputdir,'Freq');
    cfg.savefile = fullfile(savepath,['tfr2_',subj]);
    [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
    
    disp(['peaked timing: ',num2str(time_of_interest),' Sec'])
    disp(['peaked freq: ',num2str(freq_of_interest),' Hz']);
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    % data.sens = sens;
    
    %% Epoching
    %     switch task
    %         case 1 %'DefNam'
    %             toi = input('Eneter toi (e.g. [-0.3,0;1.5,2]): ');
    %         case 2 %'PicNam'
    %             %             toi = [-0.3,0;0.7,1];
    %             toi = [-0.2,0;0.5,0.8];
    %     end
%     toi = [-0.2,0;1,1.5];
    toi = [-0.3,0;0.7,1.5];
    ep_data = vy_epoch(cln_data, toi);
    
    %- Appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Data-covarinace estimation
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');
    
    %% Grand Mean
    a_data = vy_ave(cln_data);
    savepath = fullfile(outputdir,'Timelock');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    cfg = [];
    cfg.savefile = fullfile(savepath,['gmean_',subj,'.mat']);
    cfg.saveflag = 1;
    cfg.lay  = lay;
    vy_ave_plot(cfg, a_data);
    
    %% Volumetric-based analysis
    mridir = fullfile(datadir,subj,'brainstorm_db/anat');
    d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
    clear fid
    if ~isempty(d)
        sMRI1 = d.name;
        load(sMRI1);
        fid.SCS = SCS;
        fid.NCS = NCS;
        mripfile = fullfile(mridir,'T1.nii');
        outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
        if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
        
        cfg = [];
        cfg.megdata = t_data.app;
        cfg.mripfile = mripfile;
        cfg.hsfile = datafile; % headshape;
        cfg.fid = fid;
        cfg.outputmridir = outputmridir;
        cfg.subj = subj;
        cfg.plotflag = 2;
        [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg);
    end
    
    %% Choosing mesh
    switch meshgrid
        case 1
            meshtag = 'lowres';
            load temp_grid % low-res
            individual_grid = individual_grid_10mm;
        case 2
            meshtag = 'highres';
            load temp_grid_8mm % high-res
            individual_grid = individual_grid_8mm;
    end
    
    %% Anatomoy check!
    anatomy_check_flag = 1;
    saveflag = 1;
    if anatomy_check_flag == 1
        vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
    end
    %     close all
    
    %%
    outputdir = fullfile(outputdir,mtag);
    %%
    switch method
        case 1
            %%
            vy_source_lcmv_light
            %                             vy_source_lcmv
        case 2
            %%
            cfg = [];
            cfg.p.ft_old = ft_old;
            cfg.p.ft_path = ft_path;
            cfg.p.cd_org = cd_org;
            cfg.p.hcp_path = hcp_path;
            cfg.grid = individual_grid;
            cfg.headmodel = individual_headmodel;
            cfg.subj = subj;
            cfg.sens = sens;
            cfg.outputdir = outputdir;
            cfg.template_grid = template_grid;
            cfg.template_mri = template_mri;
            vy_network_light1(cfg,t_data) % conn-network analysis
        case 3
            %%
            cfg = [];
            cfg.grid = individual_grid;
            cfg.headmodel = individual_headmodel;
            cfg.sens = sens;
            cfg.outputdir = outputdir;
            cfg.template_grid = template_grid;
            cfg.template_mri = template_mri;
            vy_source_dics(cfg,ep_data);
            
        case 4
            %%
            outputdir = fullfile(outdir,'ft_process',yttag, subj, tag);
            outputdir1 = fullfile(outputdir, 'spm_source');
            if exist(outputdir1, 'file') == 0
                mkdir(outputdir1);   % create a directory
            end
            
            cfg = [];
            cfg.toilim = [-0.4 2];
            eint_data = ft_redefinetrial(cfg, cln_data);
            
            if exist(mripfile, 'file') == 2
                cd(outputdir1);
                
                cfg = [];
                cfg.p.spm = spm_path;
                cfg.p.hcp_path = hcp_path;
                cfg.p.ft_path = ft_path;
                cfg.p.spmbf = spmbf_path;
                cfg.p.cd_org = cd_org;
                cfg.datafile = datafile;
                cfg.eint_data = eint_data;
                cfg.mripfile = mripfile;
                cfg.subj = subj;
                vy_forward_spm_meg(cfg);
                
                restoredefaultpath
                addpath(genpath(ft_path));
                addpath(genpath(hcp_path));
                addpath(genpath([cd_org,'/functions']));
                addpath(genpath([cd_org,'/Data_file']));
                cd(outputdir)
            end
    end
    
    %%
    pause
    close all
    disp([datafile,' ,was completed'])
end





