clear; clc, close('all'); warning off

%% Initial settings
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));
vy_init

outputdir1 = '/home/vyoussof/Desktop/Vahab/Processed_data'; % saving directory
datadir = '/home/vyoussof/Desktop/Vahab';

epoch_type = 'STI101';

%%
subj = 'FD'; run = '1';
mridir = fullfile(datadir,subj);

% % - Auditory definition naming
% task = 'DefNam';
% datafolder = fullfile(mridir,'Run09_DFNM_vertical');
% datafile = fullfile(datafolder,'Run09_DFNM_vertical_raw_tsss.fif');

% - Visual picture naming
task = 'PicNam';
datafolder = fullfile(mridir,'Run08_PN_vertical');
datafile = fullfile(datafolder,'Run08_PN_vertical_raw_tsss.fif');

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
        Evnt_IDs = [2,3]; % 3: images, 2: scrambled images
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
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outputdir,['beta_f_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        f_data = vy_preprocess_beta(datafile,Evnt_IDs,epoch_type);
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
    %% ICA cleaning
    savepath = ['beta_ica_',subj,'_',run,'.mat'];
    saveflag = 1;
    ic_data = vy_ica_cleaning(f_data, lay, savepath, saveflag);
    
    %%
    disp('rej bad /channels/trials ...');
    savepath = fullfile(outputdir,['beta_r_',subj,'_',run,'.mat']);
    saveflag = 2;
    cln_data = vy_rejectvisual(ic_data, lay, savepath, saveflag);
    
    %% Data interpolation
    %     load('neuromag306mag_neighb');
    %     [data_clean, badchans] = vy_interpolate_meg(data_fix, f_data, neighbours, 0);
    
    %% data inspection
    %     cfg = [];
    %     cfg.viewmode = 'vertical';
    %     cfg.continuous = 'yes';
    %     ft_databrowser(cfg,data_clean);
    
    %% freq analysis (fft)
    savepath = fullfile(outputdir,['beta_fft_',subj,'_',run]);
    vy_fft(cln_data, [2 40],1,savepath,1);
    % vy_fft(data_clean, [2 40],1,[],0);
    
    % pausue
    %% freq analysis (tfr)
    savepath = fullfile(outputdir,['beta_tfr_ica',subj,'_',run]);
    vy_tfr(cln_data,lay,savepath);
    %     vy_tfr(data_clean,lay,[]);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    % data.sens = sens;
    
    %% Epoching
    %     toi = [-0.3,0;0.4,0.8];
    toi = [0,0.3;0.5,0.8];
    %     toi = [0, 0.3;0.5,0.8];
    ep_data = vy_epoch(cln_data, toi);
    
    %% appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Data-covarinace estimation
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');
    
    %% Grand Mean
    a_data = vy_ave(cln_data);
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
            outputdir_lcmv = fullfile(outputdir,'lcmv');
            vy_source_lcmv_light
            
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
    
    %%
    %     cfg = [];
    %     cfg.chantype    = 'megmag';
    %     cfg.preproc.baselinewindow = [-inf 0];
    %     cfg.trials = find(f_data.trialinfo==3); % images
    %     data_img = ft_redefinetrial(cfg, f_data);
    %
    %     cfg.trials = find(f_data.trialinfo==2); % Scrambled images
    %     data_scimg = ft_redefinetrial(cfg, f_data);
    
    %     %% ICA cleaning
    %     savepath = ['ica_imag',subj,'_',run,'.mat'];
    %     data_img_fix = vy_ica_cleaning(data_img, lay, savepath, 2);
    %
    %     savepath = ['ica_scram',subj,'_',run,'.mat'];
    %     data_scram_fix = vy_ica_cleaning(data_scimg, lay, savepath, 2);
    %
    %     %%
    %     savepath = ['ica_imag',subj,'_',run,'.mat'];
    %     data_img_fix = vy_ica_cleaning(data_img, lay, savepath, 2);
    %
    %     savepath = ['ica_scram',subj,'_',run,'.mat'];
    %     data_scram_fix = vy_ica_cleaning(data_scimg, lay, savepath, 2);
    %     %%
    %     data_img_fix1 = vy_rejectvisual(data_img_fix, lay);
    %     data_scram_fix1 = vy_rejectvisual(data_scram_fix, lay);
    %
    % %
    %     %% Data interpolation
    %     %     load('neuromag306mag_neighb');
    %     %     [data_clean, badchans] = vy_interpolate_meg(data_fix, f_data, neighbours, 0);
    % %     data_clean = data_fix;
    %
    %     %% data inspection
    %     %     cfg = [];
    %     %     cfg.viewmode = 'vertical';
    %     %     cfg.continuous = 'yes';
    %     %     ft_databrowser(cfg,data_clean);
    %
    %     %% freq analysis (fft)
    %     %     savepath = fullfile(outputdir,['fft_',subj,'_',run]);
    %     vy_fft(data_img_fix1, [2 40],1,[],0);
    %     vy_fft(data_scram_fix1, [2 40],1,[],0);
    %
    %     % pausue
    %     %% freq analysis (tfr)
    %     %     savepath = fullfile(outputdir,['tfr_ica',subj,'_',run]);
    %     %     vy_tfr(data_clean,lay,savepath);
    %     vy_tfr(data_img_fix1,lay,[]);
    %     vy_tfr(data_scram_fix1,lay,[]);
    %
    %     %% elec/grad
    %     sens = ft_read_sens(datafile);
    %     sens = ft_convert_units(sens,'mm');
    %     % data.sens = sens;
    %
    %     %% Epoching
    %     toi = [0, 0.3; 0,0.8];
    %     cfg = [];
    %     cfg.toilim = toi(2,:);
    %     ep_img = ft_redefinetrial(cfg, data_img_fix1);
    %     ep_scrm = ft_redefinetrial(cfg, data_scram_fix1);
    %
    %     %%
    %     cfg                  = [];
    %     cfg.covariance       = 'yes';
    %     cfg.covariancewindow = 'all';
    %     cfg.preproc.demean   = 'yes';    % enable demean to remove mean value from each single trial
    %     cfg.keeptrials       = 'yes';
    %     t_img                = ft_timelockanalysis(cfg, ep_img);
    %     t_scrm               = ft_timelockanalysis(cfg, ep_scrm);
    %
    %     %% Grand Mean
    %     a_data = vy_ave(t_img);
    %     vy_ave_plot(a_data,lay,[])
    %
    %     a_data = vy_ave(t_scrm);
    %     vy_ave_plot(a_data,lay,[])
    %
    %     %% mri anatomy
    %     mripath = fullfile(mridir,'T1.nii');
    %     outputmridir = fullfile(outputdir1,subj,'anat'); % output dir
    %
    %     if exist(outputmridir, 'file') == 0
    %         mkdir(outputmridir);   %create the directory
    %     end
    %     hsfile = datafile; % headshape
    %
    %     if exist(mripath, 'file') == 2
    %
    %         [mri_realigned,individual_seg, individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj);
    %         if exist('mri_realigned', 'var') == 1
    %
    %             %- high res grid
    %             load temp_grid_8mm % from, vy_warping()
    %
    %             %% Source modelling (warped with template)
    %             cfg                 = [];
    %             cfg.grid.warpmni    = 'yes';
    %             cfg.grid.nonlinear  = 'yes';
    %             cfg.grid.template   = template_grid;
    %             cfg.mri             = mri_realigned;
    %             cfg.grid.unit       = 'mm';
    %             individual_grid     = ft_prepare_sourcemodel(cfg);
    %
    %             %% anatomoy check!
    %             %     vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir)
    %             %             vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,[])
    %
    %             %% lcmv source analysis - whole-brain (warped with template)
    %             %             outputdir_lcmv = fullfile(outputdir,'lcmv');
    %             %             vy_source_lcmv
    %
    %             %% lcmv source analysis - whole-brain (warped with template)
    %             %             outputdir_lcmvstat = fullfile(outputdir,'lcmvstat');
    %             %             vy_source_lcmv_stat
    %
    %             %% dics (18-24Hz) source analysis - whole-brain (warped with template)
    %             %             outputdir_dics = fullfile(outputdir,'dics');
    %             %             vy_source_dics
    %
    %             %% conn & network analysis
    %             t_data.bsl = t_img;
    %             t_data.pst = t_scrm;
    %
    %             outputdir_net = fullfile(outputdir,'network_test');
    %             vy_network_light1
    %
    %             %%
    %             clc
    %             close all
    %             disp([datafile,' ,was completed'])
    %         else
    %             disp([datafile,' ,was not completed'])
    %         end
    %     end
    
end




