clear; clc, close('all'); warning off

%% initial settings
restoredefaultpath
addpath(genpath('.\functions'));
addpath .\Data_file;

ft_path = 'F:\My Matlab\My codes\My GitHub\fieldtrip_041718\fieldtrip-master';
addpath(ft_path);
ft_defaults % this loads the rest of the defaults;

hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
addpath(genpath(hcp_path));

atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal\ROI_MNI_V4.nii'));

load temp_grid_8mm % from, vy_warping()
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

DestDirectory = 'H:\VNS'; % saving directory

ft_old = 'F:\My Matlab\My codes\My GitHub\fieldtrip';

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

%%
for i = 20:20% size(datafolder,1)
    
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
    savepath = fullfile(outputdir,['r_',subj,'_',run,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        disp('preprocessing was completed');
        % savepath = fullfile(outputdir,['f_',subj,'_',run,'.mat']);
        % save(savepath, 'f_data', '-v7.3');
        
        disp('Artifact rej ...');
        cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.layout   = lay;   % this allows for plotting individual trials
        cfg.latency = [0,900];
        r_data   = ft_rejectvisual(cfg, f_data);
        save(savepath, 'r_data', '-v7.3');
        textfile_rej = fullfile(outputdir,'rej');
        rej = r_data.cfg.artfctdef.summary.artifact;
        hcp_write_ascii(textfile_rej, 'rej');
    end
    %% ica (step 1)
    comp = vy_ica(r_data,lay);
    title(subj)
    
    %- rej component
    cfg = [];
    cfg.updatesens = 'no';
    %     cfg.component = comp.label(51:end);
    bic = input('Select bad ICs:');
    cfg.component = comp.label(bic);
    data_fix = ft_rejectcomponent(cfg, comp, r_data);
    
    %% ica (step 2)
    comp = vy_ica(data_fix,lay);
    title(subj);
    
    %- rej component
    cfg = [];
    cfg.updatesens = 'no';
    %     cfg.component = comp.label(51:end);
    bic = input('Select bad ICs:');
    cfg.component = comp.label(bic);
    data_fix = ft_rejectcomponent(cfg, comp, data_fix);
    
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
    % vy_tfr(r_data,lay);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    % data.sens = sens;
    
    %% Epoching
    toi = [-0.4, -0.1;0.4,0.7];
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
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    outputmridir = fullfile(DestDirectory,'Processed','MEG','anat', subj); % output dir
    
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = fullfile(datafolder1,'hs_file'); % headshape
    
    [mri_realigned,individual_seg,~,individual_headmodel,headshape] = vy_mri(mripath,hsfile,outputmridir,subj);
    
    if exist('mri_realigned', 'var') == 1
        
        %% mesh-grid
        % - low res mesh grid
        load temp_grid
        template_grid = ft_convert_units(template_grid, 'mm');
        
        %- high res grid
        %     load temp_grid_8mm % from, vy_warping()
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
        vy_network
%         vy_network_light
        
        %% conn & network analysis - single trial
        %         vy_network_st
        
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




