clear; clc, close('all')

cd 'F:\My Matlab\My codes\My scripts\UTHSC\meg\My scripts\ft_updated'
cd_org = cd;
addpath(genpath(fullfile(cd_org,'functions')));
vy_init(cd_org)
spm_path = spm('dir');

hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

addpath('F:\My Matlab\My codes\My GitHub\fieldtrip_041718\fieldtrip-master\external\eeglab');

%%
DestDirectory = 'H:\VNS'; % saving directory
ft_old = 'F:\My Matlab\My codes\My GitHub\fieldtrip';
mridir = 'H:\VNS\MRI\Nifti\T1';
p = fullfile(DestDirectory,'MEG');

disp('1: CRM');
disp('2: VG-Auditory')
disp('3: VG-Printed')
task = input('Enter task number: ');

switch task
    case 1
        tsk = 'CRM'; % Continuous recognition memory
        name = 'CRM';
    case 2
        tsk = 'VGA'; % VerbGen-Auditory
        name = 'VerbGenAud';
    case 3
        tsk = 'VGP'; % VerbGen-Printed
        name = 'VerbGenVis';
end

d = rdir([p,'\**\',name,'\**\c,rfDC']);
% d = rdir([p,'\**\',task,'\**\config']);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
end
datafolder = datafolder';
disp(datafolder)

%% ICA cleaned data
disp('listing ICA cleaned data ...');
data_preprocess = [];
p = fullfile('H:\VNS\Preprocessed',tsk);
d = rdir([p, '\*.mat']);
for i=1:length(d)
    data_preprocess{i} = d(i).name;
end
data_preprocess = data_preprocess';
disp(data_preprocess),

%% Processing method
disp('1: Source_lcmv');
disp('2: Source_lcmv-stats');
disp('3: Source_dics')
disp('4: Network_evc')
disp('5: SPM-COH(loreta)');
method = input('Eneter the method: ');

%%
cd_org = cd;
for i = 2:size(datafolder,1)
    
    datafolder1 = datafolder{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    datafile = fullfile(datafolder1,'c,rfDC');
    Index = strfind(datafolder1, '\');
    subj = datafolder1(Index(3)+1:Index(4)-1);
    run  = datafolder1(Index(5)+1:end);
    disp(['subj:',subj])
    disp(['run:',run])
    outputdir = fullfile(DestDirectory,'Processed','MEG','ft3',tsk, subj,run);
    disp(['outputdir:',outputdir])
    if exist(outputdir, 'file') == 0
        mkdir(outputdir);   %create the directory
    end
    disp([datafolder1,' is analysing']),
    cd(outputdir)
            
    %% 4D layout
    cfg = [];
    cfg.layout = '4D248.mat';
    lay = ft_prepare_layout(cfg);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outputdir,['f_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
%         f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        f_data = vy_preprocess5(datafile,task);
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
     %% reading cleaned data (with ICA)
%     for k =1:length(data_preprocess)
%         datafolder2 = data_preprocess{k}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
%         %     datafile = fullfile(datafolder1,tsk);
%         Index = strfind(datafolder2, tsk);
%         subj1 = datafolder2(Index(1)+4:Index(1)+8);
%         switch task
%             case 1
%                 run1  = datafolder2(Index(2)+6);
%             case 2
%                 run1  = datafolder2(Index+21);
%             case 3
%                 run1  = datafolder2(Index+21);
%         end
%         if strcmp(subj,subj1) && run == run1
%             disp(['subj:',subj]);
%             disp(['run:',run]);
%             disp('loading ... ');
%             disp(data_preprocess{k});
%             load(data_preprocess{k});
%             disp('matched with ...')
%             disp(datafile);
%             break,
%         end
%     end
%     if (strcmp(subj,subj1) && run == run1) == 0
%         error(['no match was found for ',subj, ' run: ', run]);
%     end
%     
%     %% Appending cleaned data, another filtering .1-40 Hz
%     savepath1 = ['ica_',subj,'_',run,'.mat'];
%     if exist(savepath1, 'file') == 2
%         load(savepath1);
%     else
%         [icdata, av_data] = vy_cleandataICA(datafile,clean_data,task);
%         save(savepath1, 'icdata', '-v7.3');
%     end
%     
%     %% rejecting bad channels and trials, as informed by ICA
%     %- chan
%     badchannels1 = setdiff(f_data.label,icdata.label);
%     disp('bad channels:');
%     disp(badchannels1)
%     badchannels2 = strcat('-',badchannels1);
%     badchannels2 = [badchannels2;'all'];
%     
%     % -trials
%     [I,badtrials] = setdiff(f_data.sampleinfo(:,1),icdata.sampleinfo(:,1));
%     disp('Bad trials:');
%     disp(badtrials)
%     [idx,loc] = ismember(f_data.sampleinfo(:,1),icdata.sampleinfo(:,1));
%     goodtrl = loc(idx);
%     
% %     %- data inspection
% %     cfg = [];
% %     cfg.viewmode = 'vertical';
% %     cfg.continuous = 'no';
% %     cfg.channel = badchannels1;
% %     ft_databrowser(cfg,f_data);
%    
%     cfg = [];
%     cfg.channel = badchannels2;
%     cfg.trials = goodtrl;
%     f_data1 = ft_selectdata(cfg, f_data);
    
    %% Rejecting bad epoches
    disp('Artifact rejection ...');
    savepath2 = ['r_',subj,'_',run,'.mat'];
    if exist(savepath2, 'file') == 2
        load(savepath2)
    else
        disp('Artifact rej ...');
        cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.layout   = lay;   % this allows for plotting individual trials
        cfg.latency = [0,900];
        r_data   = ft_rejectvisual(cfg, f_data);
        save(savepath2, 'r_data', '-v7.3');
    end
    
    %% Downsampling
    cfg = [];
    cfg.resamplefs = 300;
    cfg.detrend = 'no';
    r_data = ft_resampledata(cfg,r_data);
    
    %% Interpolating data
    disp('Interpolating data ...');
    savepath4 = ['int_ica_',subj,'_',run,'.mat'];
    if exist(savepath4, 'file') == 2
        load(savepath4)
    else
        [int_data, badchans] = vy_interpolate_meg(r_data, f_data);
        save(savepath4, 'int_data', '-v7.3');
        disp('Interpolated chans:')
        disp(badchans);
    end
    
    %%
%     data_clean = int_data;
    data_clean = r_data;
    
    %% Grand Mean
    a_data = vy_ave(data_clean);
    savepath = fullfile(outputdir,['gmean_',subj,'_',run]);
    t_max = vy_ave_plot(a_data,lay,savepath);
%     vy_ave_plot(a_data.all,lay,[])
    
    %% freq analysis (fft)
%     savepath = fullfile(outputdir,['fft_ica',subj,'_',run]);
%     vy_fft(data_clean, [2 40],1,savepath,1);
    %     vy_fft(data_clean, [2 40],1,[],0);
    
    %% freq analysis (tfr)
    %     savepath = fullfile(outputdir,['tft_ica',subj,'_',run]);
    %     vy_tfr(data_clean,lay,savepath);
    
    %% Epoching
    toi = [-0.4, -0.1;0.4,0.7];
    data_clean.grad = sens;
    ep_data = vy_epoch(data_clean,toi);
    
    %% appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Timelock
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');

    %% mri anatomy   
    switch subj
        case 'TQuality1'
            subj = 'C-117';
        case 'TQuality2'
            subj = 'C-120';
    end
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    outputmridir = fullfile(DestDirectory,'Processed','MEG','anat', subj); % output dir
    vy_init(cd_org)
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = fullfile(datafolder1,'hs_file'); % headshape
    [mri_realigned,~,~,individual_headmodel,~] = vy_mri(mripath,hsfile,outputmridir,subj);
    
    [individual_headmodel, sens] = ft_prepare_vol_sens(individual_headmodel, sens, 'channel', ep_data.all.label);
    vy_init(cd_org)
    
    %% mesh-grid
    % - low res mesh grid
    load temp_grid
    template_grid = ft_convert_units(template_grid, 'mm');
    
    %- high res grid
%     load temp_grid_8mm % from, vy_warping()
    
    %- Warpiing to a template with 8 mm res
    cfg                 = [];
    cfg.grid.warpmni    = 'yes';
    cfg.grid.nonlinear  = 'yes';
    cfg.grid.template   = template_grid;
    cfg.mri             = mri_realigned;
    cfg.grid.unit       = 'mm';
    individual_grid     = ft_prepare_sourcemodel(cfg);
    
    %% anatomy inspections ( headmodel, mesh, alignment, ...)
    %  vy_mri_inspection(volp,individual_grid,headshape, mri_realigned,sens, outputmridir)
    
    %% perform whole-brain source reconstruction
    switch method
        case 1
            mtd = 'lcmv';
            vy_source_lcmv
        case 2
            mtd = 'lcmv-stas';
            vy_source_lcmv_stat
        case 3
            mtd = 'dics';
            disp('wavelet analysis ...');
            savepath5 = ['wavelet_',subj,'_',run,'.mat'];
            if exist(savepath5, 'file') == 2
                load(savepath5)
            else
                savepath = fullfile(outputdir,['wavelet_',subj,'_',run]);
                foi = 1:1:30;
                w_data = vy_wavelet(data_clean, lay, savepath,foi);
            end           
            vy_source_dics2
            vy_source_dics
        case 4            
            mtd = 'net';
            vy_network_light           
        case 5
            cfg = [];
            cfg.toilim = [-0.4 1];
            eint_data = ft_redefinetrial(cfg, data_clean);
            
            mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
            name = [subj,'_',run];
            
            outputdir1 = fullfile(outputdir, 'spm_source');
            if exist(outputdir1, 'file') == 0
                mkdir(outputdir1);   % create a directory
            end
            cd(outputdir1);
            vy_init(cd_org)
            vy_forward_spm_meg(datafile,eint_data,mripath,subj);
            cd(outputdir)
    end
end
% end
 

% %% look at the analysis history
% prefix = sprintf('/tmp/Sub%02d', subj);
% cfg           = [];
% cfg.filename  = [prefix '_avg_Faces_vs_Scrambled.html'];
% ft_analysispipeline(cfg, avg_Faces_vs_Scrambled);



