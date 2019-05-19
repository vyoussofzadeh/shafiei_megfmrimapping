clear; clc, close('all')

cd 'F:\My Matlab\My codes\My scripts\UTHSC\meg\My scripts\ft_updated'
cd_org = cd;
addpath(genpath(fullfile(cd_org,'functions')));
vy_init(cd_org)


hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
load temp_grid_8mm % from, vy_warping()
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

%%
DestDirectory = 'H:\VNS'; % saving directory
ft_old = 'F:\My Matlab\My codes\My GitHub\fieldtrip';
mridir = 'H:\VNS\MRI\Nifti\T1';

disp('1: CRM');
disp('2: VG-Auditory')
disp('3: VG-Printed')
task = input('Enter task number: ');
switch task
    case 1
        tsk = 'CRM'; % Continuous recognition memory
        if exist('CRM_datafolders.mat','file')== 2
            load CRM_datafolders
        else
            datafolder = spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
            save('CRM_datafolders2','datafolder');
        end
        disp(datafolder)
    case 2
        tsk = 'VGA'; % VerbGen-Auditory
        if exist('VGA_datafolders.mat','file')== 2
            load VGA_datafolders
        else
            datafolder = spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\VGA\1
            save('VGA_datafolders','datafolder');
        end
        disp(datafolder)
    case 3
        tsk = 'VGP'; % VerbGen-Printed
        if exist('VGP_datafolders.mat','file')== 2
            load VGP_datafolders
        else
            datafolder = spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\VGP\1
            save('VGP_datafolders','datafolder');
        end
        disp(datafolder)
end

cd_org = cd;
%%
for i=1:size(datafolder,1)
    
    % Data select
    datafolder1 = datafolder(i,:); % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    datafile = fullfile(datafolder1,'c,rfDC');
    Index = strfind(datafolder1, '\');
    subj = datafolder1(Index(3)+1:Index(4)-1);
    run  = datafolder1(Index(5)+1:end);
    disp(['subj:',subj])
    disp(['run:',run])
    outputdir = fullfile(DestDirectory,'Processed','MEG','ft',tsk, subj,run);
    disp(['outputdir:',outputdir])
    if exist(outputdir, 'file') == 0
        mkdir(outputdir);   %create the directory
    end
    disp([datafolder1,' is analysing']),
    cd(outputdir)
    
    %% trial check
    %     vy_trl_check(datafile,task)
    
    %% reading events
    event = ft_read_event(datafile);
    % Evnt_IDs = {event(5).value};
    switch tsk
        case 'CRM'
            epoch_type = {'TRIGGER'};
            %             Evnt_IDs = {event(4).value};
            Evnt_IDs = {128};
        case 'VGA'
            epoch_type = {'TRIGGER'};
            Evnt_IDs = {64};
            
        case 'VGP'
            epoch_type = {'TRIGGER'};
            %             Evnt_IDs = {event(5).value};
            Evnt_IDs = {2052};
    end
    disp([Evnt_IDs,epoch_type])
    
    %% hcp check
    % currentFolder = pwd;
    % qadir = fullfile(outputdir,'qa');
    % if exist(qadir, 'file') == 0
    %     mkdir(qadir);   %create the directory
    % end
    % hcp_datacheck_vy(datafile,qadir)
    % % vy_datacheck(datafile)
    % cd(currentFolder);
    
    %% 4D layout
    cfg = [];
    % cfg.layout = '4D248.lay';
    cfg.layout = '4D248.mat';
    lay = ft_prepare_layout(cfg);
    % ft_layoutplot(cfg);
    % pause
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outputdir,['r_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        %         f_data = vy_preprocess3(datafile,Evnt_IDs,epoch_type,task);
        disp('preprocessing was completed');
        % savepath = fullfile(outputdir,['f_',subj,'_',run,'.mat']);
        % save(savepath, 'f_data', '-v7.3');
        
        disp('Artifact rej ...');
        cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.layout   = lay;   % this allows for plotting individual trials
        r_data   = ft_rejectvisual(cfg, f_data);
        save(savepath, 'r_data', '-v7.3');
        textfile_rej = fullfile(outputdir,'rej');
        rej = r_data.cfg.artfctdef.summary.artifact;
        hcp_write_ascii(textfile_rej, 'rej');
    end
    
    %%
    cfg = [];
    cfg.method = 'pca';
    cfg.updatesens = 'no';
    comp = ft_componentanalysis(cfg, r_data);
    
    %see the components and find the artifact
    cfgb=[];
    cfgb.layout = lay;
    cfgb.channel = {comp.label{1:5}};
    cfgbo=ft_databrowser(cfgb,comp);
    
    cfg = [];
    cfg.updatesens = 'no';
    cfg.component = comp.label(51:end);
    data_fix = ft_rejectcomponent(cfg, comp);
    r_data2 = data_fix;
    %         ft_analysispipeline([], r_data);
    
    %% ica
    % i_data = vy_ica(data,lay);
    
    %% data inspection
    % cfg = [];
    % cfg.viewmode = 'vertical';
    % cfg.continuous = 'yes';
    % ft_databrowser(cfg,f_data);
    
    %% freq analysis (fft)
    %     savepath = fullfile(outputdir,['fft_',subj,'_',run]);
    %     vy_fft(r_data, [2 40],1,savepath,1);
    %     vy_fft(r_data, [2 40],1,[],0);
    
    % pausue
    %% freq analysis (tfr)
    %     savepath = fullfile(outputdir,['tft_',subj,'_',run]);
    %     vy_tfr(r_data,lay,savepath);
    
    %% elec/grad
    %     sens = ft_read_sens(datafile);
    % data.sens = sens;
    %% spm analysis
    cfg                         = [];
    cfg.dataset                 = datafile;
    cfg.channel = {'MEG'};
    raw_data = ft_preprocessing(cfg);
    
    [int_data, badchans] = vy_interpolate_meg(r_data2, raw_data);
    
    cfg = [];
    cfg.toilim = [-0.4 1];
    eint_data = ft_redefinetrial(cfg, int_data);
    
    outputdir1 = fullfile(outputdir, 'spm_source');
    if exist(outputdir1, 'file') == 0
        mkdir(outputdir1);   %create the directory
    end
    cd(outputdir1);
    
    %     sMRI = fullfile(mridir,[subj,'_T1.nii']);
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    if exist(mripath, 'file') ~= 2
        sMRI = 'F:\My Matlab\SPM\spm12_4\spm12\spm12\canonical\single_subj_T1.nii';
    end
    
    name = [subj,'_',run];
    
    %     vy_forward_spm_meg(datafile,eint_data,mripath,subj);
    vy_source_spm_meg2(datafile,eint_data,mripath,subj)
    vy_init(cd_org)
    
    %% Epoching
    ep_data = vy_epoch(int_data);
    
    %% appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Timelock
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');
    %
    %     %% Grand Mean
    %     a_data = vy_ave(ep_data);
    %     savepath = fullfile(outputdir,['gmean_',subj,'_',run]);
    %     vy_ave_plot(a_data.all,lay,savepath)
    %     %     vy_ave_plot(a_data.all,lay,[])
    %
    %% mri anatomy
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    outputmridir = fullfile(DestDirectory,'Processed','MEG','anat', subj); % output dir
    vy_init(cd_org)
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = fullfile(datafolder1,'hs_file'); % headshape
    [mri_realigned,individual_seg,~,individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj);
    
    [volp, sens] = ft_prepare_vol_sens(individual_headmodel, sens, 'channel', ep_data.all.label);
    vy_init(cd_org)
    
    %% low res grid
    load temp_grid
    template_grid = ft_convert_units(template_grid, 'mm');
    
    %% high res grid
    % load temp_grid_8mm % from, vy_warping()
    
    %%
    % warpiing with new template res
    cfg                 = [];
    cfg.grid.warpmni    = 'yes';
    cfg.grid.nonlinear  = 'yes';
    cfg.grid.template   = template_grid;
    cfg.mri             = mri_realigned;
    cfg.grid.unit       = 'mm';
    individual_grid     = ft_prepare_sourcemodel(cfg);
    
    %% mr inspection
    % vy_mri_inspection(volp,individual_grid,headshape, mri_realigned,sens, outputmridir)
    % vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,[])
    
    %% lcmv source analysis - whole-brain (warped with template)
    vy_source_lcmv
    
    %% lcmv source analysis - whole-brain (warped with template)
    %         vy_source_lcmv_stat
    
    %% dics (18-24Hz) source analysis - whole-brain (warped with template)
    %     vy_source_dics
    
    %% conn & network analysis
    vy_network
    %         clear s_data2
    
    %% conn & network analysis - single trial
    %         vy_network_st
    
    %% conn & network analysis - ImCoh (spectral)
    %         vy_network2
    
    %         clc
    %         close all
    %         disp([datafile,' ,was completed'])
    %     else
    disp([datafile,' ,was not completed'])
end

% end

%% look at the analysis history
prefix = sprintf('/tmp/Sub%02d', subj);
cfg           = [];
cfg.filename  = [prefix '_avg_Faces_vs_Scrambled.html'];
ft_analysispipeline(cfg, avg_Faces_vs_Scrambled);



