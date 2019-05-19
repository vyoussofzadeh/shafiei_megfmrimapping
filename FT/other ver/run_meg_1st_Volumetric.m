clear; clc, close('all'); warning off

%% Initial settings
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));
vy_init

outputdir1 = '/home/vyoussof/Desktop/Vahab/Processed_data'; % saving directory
datadir = '/home/vyoussof/Desktop/Vahab';

%%
subj = 'FD'; run = '1'; 
mridir = fullfile(datadir,subj);

% % - Auditory definition naming 
task = 'DefNam';
datafolder = fullfile(mridir,'Run09_DFNM_vertical');
datafile = fullfile(datafolder,'Run09_DFNM_vertical_raw_tsss.fif');

% - Visual picture naming
% task = 'PicNam';
% datafolder = fullfile(mridir,'Run08_PN_vertical');
% datafile = fullfile(datafolder,'Run08_PN_vertical_raw_tsss.fif');

%%
% subj = 'PB'; run = '1';
% mridir = fullfile(datadir,subj);
% 
% % - Auditory definition naming
% task = 'DefNam';
% datafolder = fullfile(mridir,'Run08_DFNAM_upright');
% datafile = fullfile(datafolder,'Run08_DFNAM_upright_raw_tsss.fif');
% 
% % - Visual picture naming, 2(=60) is scrambled images, 3(=120) is images. 
% task = 'PicNam';
% datafolder = fullfile(mridir,'Run06_PN_1_2_upright');
% datafile = fullfile(datafolder,'Run06_PN_1_2_upright_raw_tsss.fif');

switch task
    case 'DefNam'
        Evnt_IDs = 1; % questions
    case 'PicNam'
        Evnt_IDs = 3; % questions
end

%%
outputdir = fullfile(outputdir1,subj, task);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end

%%
for i = 1:size(datafolder,1)
    
    disp([datafolder,' is analysing']),
    cd(outputdir)
    
    %% 4D layout
    cfg = [];
    % cfg.layout = '4D248.lay';
    cfg.layout = 'neuromag306mag.lay';
    lay = ft_prepare_layout(cfg);
    % ft_layoutplot(cfg);
    % pause
    
    %%
%     event = ft_read_event(datafile);
    
    %%
    epoch_type = 'STI101';
    f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
    
    %%
%     cfg = [];
%     cfg.method = 'pca';
%     cfg.updatesens = 'no';
%     cfg.chantype    = 'megmag';
%     comp = ft_componentanalysis(cfg, f_data);
%     
%     cfg = [];
%     cfg.updatesens = 'no';
%     cfg.component = comp.label(51:end);
%     data_fix = ft_rejectcomponent(cfg, comp);
    
%     data_fix = f_data;
%     
%     %%
%     cfg = [];
%     cfg.chantype    = 'megmag';
%     cfg.preproc.baselinewindow = [-inf 0];
%     cfg.trials = find(data_fix.trialinfo==1); % Questions
%     data_ques = ft_redefinetrial(cfg, data_fix);
    
%     cfg.trials = find(data_fix.trialinfo==2); % Scrambled Audio
%     data_audsc = ft_redefinetrial(cfg, data_fix);
    
    %% Filteting, Event reading, Artifact rejecrtion
%     disp('preprocessing ...');
%     savepath = fullfile(outputdir,['f_',subj,'_',run,'.mat']);
%     if exist(savepath, 'file') == 2
%         load(savepath)
%     else
%         %         f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
%         f_data = vy_preprocess_meg(datafile,4);
%         disp('preprocessing was completed');
%         save(savepath, 'f_data', '-v7.3');
%     end
%     
    %% ICA cleaning
    savepath = ['ica_',subj,'_',run,'.mat'];
    
    satis = 0;
    disp('ica cleaning ...');
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.latency = [-400,900];
        cfg.layout   = lay;   % this allows for plotting individual trials
        r_data   = ft_rejectvisual(cfg, f_data);
        
        comp = vy_ica(r_data,lay);
        title(savepath)
        %rej component
        cfg = [];
        cfg.updatesens = 'no';
        bic = input(['Select bad ICs for ',subj,':']);
        cfg.component = comp.label(bic);
        data_fix = ft_rejectcomponent(cfg, comp, r_data);
        while (satis == 0)&&(isempty(bic) == 0)
            close all
            comp = vy_ica(data_fix,lay);
            satis = input('Statisfied with ICA, yes = 1, not yet = 0:');
            if satis == 1, close('all'), break, end
            bic = input(['Select bad ICs for ',subj,':']);
            cfg.component = comp.label(bic);
            data_fix = ft_rejectcomponent(cfg, comp, data_fix);
        end
        close all
        save(savepath, 'data_fix', '-v7.3');
    end
    
    %%
    disp('rej bad /channels/trials ...');
    savepath = fullfile(outputdir,['r1_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.latency = [-400,900];
        cfg.layout   = lay;   % this allows for plotting individual trials
        data_fix   = ft_rejectvisual(cfg, data_fix);
%         save(savepath, 'data_fix', '-v7.3');
    end
    
    %% Data interpolation
%     load('neuromag306mag_neighb');
%     [data_clean, badchans] = vy_interpolate_meg(data_fix, f_data, neighbours, 0);
    data_clean = data_fix;
    
    %% data inspection
    %     cfg = [];
    %     cfg.viewmode = 'vertical';
    %     cfg.continuous = 'yes';
    %     ft_databrowser(cfg,data_clean);
    
    %% freq analysis (fft)
    savepath = fullfile(outputdir,['fft_',subj,'_',run]);
    vy_fft(data_clean, [2 40],1,savepath,1);
    % vy_fft(data_clean, [2 40],1,[],0);
    
    % pausue
    %% freq analysis (tfr)
    savepath = fullfile(outputdir,['tfr_ica',subj,'_',run]);
    vy_tfr(data_clean,lay,savepath);
    %     vy_tfr(data_clean,lay,[]);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    % data.sens = sens;
    
    %% Epoching
%     toi = [-0.3,0;0.4,0.8];
    toi = [-0.3,0;1.0,2];
    %     toi = [0, 0.3;0.5,0.8];
    ep_data = vy_epoch(data_clean, toi);
    
    %% appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Data-covarinace estimation
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');
    
    %% Grand Mean
    a_data = vy_ave(data_clean);
    savepath = fullfile(outputdir,['gmean_',subj,'_',run]);
    vy_ave_plot(a_data,lay,savepath)
    
    %% mri anatomy
    mripath = fullfile(mridir,'T1.nii');
    outputmridir = fullfile(outputdir1,subj,'anat'); % output dir
    
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = datafile; % headshape
    
    if exist(mripath, 'file') == 2
        
        [mri_realigned,individual_seg, individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj);
        if exist('mri_realigned', 'var') == 1
            
            %- high res grid
%             load temp_grid_8mm % from, vy_warping()
            load temp_grid
            
            %% Source modelling (warped with template)
            cfg                 = [];
            cfg.grid.warpmni    = 'yes';
            cfg.grid.nonlinear  = 'yes';
            cfg.grid.template   = template_grid;
            cfg.mri             = mri_realigned;
            cfg.grid.unit       = 'mm';
            individual_grid     = ft_prepare_sourcemodel(cfg);
            
            %% anatomoy check!
            %     vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir)
%             vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,[])
            
            %% lcmv source analysis - whole-brain (warped with template)
%             outputdir_lcmv = fullfile(outputdir,'lcmv');
%             vy_source_lcmv
            
            %% lcmv source analysis - whole-brain (warped with template)
            %             outputdir_lcmvstat = fullfile(outputdir,'lcmvstat');
            %             vy_source_lcmv_stat
            
            %% dics (18-24Hz) source analysis - whole-brain (warped with template)
%             outputdir_dics = fullfile(outputdir,'dics');
%             vy_source_dics
            
            %% conn & network analysis
            outputdir_net = fullfile(outputdir,'network_highres');
%             vy_network
            vy_network_light1
            
            %%
            clc
            close all
            disp([datafile,' ,was completed'])
        else
            disp([datafile,' ,was not completed'])
        end
    end
    
end




