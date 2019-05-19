clear; clc, close('all'); warning off

%% initial settings
restoredefaultpath
cd '\\utdrive.uthsc.edu\babajanilab\VNS\Scripts\MEG\ft';
cd_org = cd;
addpath(genpath(fullfile(cd_org,'functions')));
vy_init(cd_org)

load temp_grid_8mm % from, vy_warping()
hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

DestDirectory = 'H:\VNS'; % saving directory

ft_old = 'F:\My Matlab\My codes\My GitHub\fieldtrip';

ft_path = 'F:\My Matlab\My codes\My GitHub\fieldtrip_packs\fieldtrip-20180809';
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal\ROI_MNI_V4.nii'));

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
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
end
datafolder = datafolder';
disp(datafolder)

%%
for i = 1:size(datafolder,1)
    
    datafolder1 = datafolder{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    datafile = fullfile(datafolder1,'c,rfDC');
    Index = strfind(datafolder1, '\');
    subj = datafolder1(Index(3)+1:Index(4)-1);
    run  = datafolder1(Index(5)+1:end);
    disp(['subj:',subj])
    disp(['run:',run])
    outputdir = fullfile(DestDirectory,'Processed','MEG','ft5',tsk, subj,run);
    disp(['outputdir:',outputdir])
    if exist(outputdir, 'file') == 0
        mkdir(outputdir);   %create the directory
    end
    disp([datafolder1,' is analysing']),
    cd(outputdir)
    
    %% check if the data have been processed!
    
    %% reading events
    event = ft_read_event(datafile);
    % Evnt_IDs = {event(5).value};
    switch tsk
        case 'CRM'
            epoch_type = event(5).type;
            Evnt_IDs = {event(4).value};
        case 'VGA'
            epoch_type = event(5).type;
            Evnt_IDs = {event(4).value};
            
        case 'VGP'
            epoch_type = event(5).type;
            Evnt_IDs = {event(5).value};
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
    
    %% Filteting, Event reading, Artifact rejecrtion
    disp('preprocessing ...');
    savepath = fullfile(outputdir,['f_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        %         f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        f_data = vy_preprocess5(datafile,task);
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end

    %% ica cleaning
    satis = 0;
    disp('ica cleaning ...');
    savepath = ['ica1_',subj,'_',run,'.mat'];
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
        bic = input(['Select bad ICs for',subj,':']);
        cfg.component = comp.label(bic);
        data_fix = ft_rejectcomponent(cfg, comp, r_data);
        while (satis == 0)&&(isempty(bic) == 0)
            close all
            comp = vy_ica(data_fix,lay);
            satis = input('Statisfied with ICA, yes = 1, not yet = 0:');
            if satis == 1, close('all'), break, end
            bic = input(['Select bad ICs for',subj,':']);
            cfg.component = comp.label(bic);
            data_fix = ft_rejectcomponent(cfg, comp, data_fix);
        end
        close all
        save(savepath, 'data_fix', '-v7.3');
    end
    %%
%     comp = vy_pca(r_data,lay);
%     bic = input(['Select bad ICs for',subj,':']);
%     cfg.component = comp.label(bic);
%     data_fix = ft_rejectcomponent(cfg, comp, r_data);
    
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
        save(savepath, 'data_fix', '-v7.3');
    end
    
    %%
    %     cfg = [];
    %     cfg.hpfilter = 'yes';
    %     cfg.hpfreq = 12;
    %     cfg.hpfiltord = 4;
    %     data_fix = ft_preprocessing(cfg,data_fix);
    
    %% data inspection
    % cfg = [];
    % cfg.viewmode = 'vertical';
    % cfg.continuous = 'yes';
    % ft_databrowser(cfg,f_data);
    
    %% freq analysis (fft)
    %     savepath = fullfile(outputdir,['fft_',subj,'_',run]);
    %     vy_fft(r_data, [2 40],1,savepath,1);
    % vy_fft(r_data, [2 40],1,[],0);
    
    % pausue
    %% freq analysis (tfr)
    % vy_tfr(data_fix,lay);
    
    %% Bandpass Filter ()
%     cfg = [];
%     cfg.bpfilter = 'yes';
%     cfg.bpfreq = [12 30];    %band-pass filter in the required range
%     data_fix = ft_preprocessing(cfg,data_fix);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    % data.sens = sens;
    
    %% Epoching
    toi = [-0.4,0;0.4,0.8];
    %     toi = [0, 0.3;0.5,0.8];
    ep_data = vy_epoch(data_fix, toi);
    
    %% appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Timelock
    t_data = vy_timelock(ep_data);
    % savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
    % save(savepath, 't_data', '-v7.3');
    
    %% Grand Mean
    a_data = vy_ave(data_fix);
    savepath = fullfile(outputdir,['gmean_',subj,'_',run]);
    vy_ave_plot(a_data,lay,savepath)
    %     vy_ave_plot(a_data,lay,[])
    
    %% mri anatomy
    switch subj
        case 'TQuality1'
            subj = 'C-117';
        case 'TQuality2'
            subj = 'C-120';
    end
    
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    outputmridir = fullfile(DestDirectory,'Processed','MEG','anat', subj); % output dir
    
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = fullfile(datafolder1,'hs_file'); % headshape
    
    if exist(mripath, 'file') == 2
        
        [mri_realigned,individual_seg,~,individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj);
        if exist('mri_realigned', 'var') == 1
            
            %% mesh-grid
            % - low res mesh grid
            load temp_grid
            template_grid = ft_convert_units(template_grid, 'mm');
            
            %- high res grid
%             load temp_grid_8mm % from, vy_warping()
            
            %% warpiing with new template res
            cfg                 = [];
            cfg.grid.warpmni    = 'yes';
            cfg.grid.nonlinear  = 'yes';
            cfg.grid.template   = template_grid;
            cfg.mri             = mri_realigned;
            cfg.grid.unit       = 'mm';
            individual_grid     = ft_prepare_sourcemodel(cfg);
            
            %% mr inspection
            %     vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir)
            % vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,[])
            
            %% lcmv source analysis - whole-brain (warped with template)
            %         vy_source_lcmv
            
            %% lcmv source analysis - whole-brain (warped with template)
            %         vy_source_lcmv_stat
            
            %% dics (18-24Hz) source analysis - whole-brain (warped with template)
            %     vy_source_dics
            
            %% conn & network analysis
%             vy_connectvitiy
            vy_network
%             vy_network_light
%             vy_network_singletrial
            
            %% conn & network analysis - single trial (stats)
            %         vy_network_st2
            
            %% conn & network analysis - ImCoh (spectral)
            %         vy_network2
            
            %%
            clc
            close all
            disp([datafile,' ,was completed'])

        else
            disp([datafile,' ,was not completed'])
        end
    end
    
end




