clear; clc, close('all'); warning off

%% initial settings
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

%%
for i = 1:size(datafolder,1)
    
    datafolder1 = datafolder{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    datafile = fullfile(datafolder1,'c,rfDC');
    Index = strfind(datafolder1, '\');
    subj = datafolder1(Index(3)+1:Index(4)-1);
    run  = datafolder1(Index(5)+1:end);
    disp(['subj:',subj])
    disp(['run:',run])
    outputdir = fullfile(DestDirectory,'Processed','MEG','ft6',tsk, subj,run);
    disp(['outputdir:',outputdir])
    if exist(outputdir, 'file') == 0
        mkdir(outputdir);   %create the directory
    end
    disp([datafolder1,' is analysing']),
    cd(outputdir)
    
    %% check if the data have been processed!
    
    %% reading events
%     event = ft_read_event(datafile);
%     % Evnt_IDs = {event(5).value};
%     switch tsk
%         case 'CRM'
%             epoch_type = event(5).type;
%             Evnt_IDs = {event(4).value};
%         case 'VGA'
%             epoch_type = event(5).type;
%             Evnt_IDs = {event(4).value};
%             
%         case 'VGP'
%             epoch_type = event(5).type;
%             Evnt_IDs = {event(5).value};
%     end
%     disp([Evnt_IDs,epoch_type])
    
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
        % cfg = [];
        cfg.metric = 'zvalue';  % use by default zvalue method
        cfg.latency = [-400,900];
        cfg.layout   = lay;   % this allows for plotting individual trials
        f_data   = ft_rejectvisual(cfg, f_data);
        
        comp = vy_ica(f_data,lay);
        title(savepath)
        %rej component
        cfg = [];
        cfg.updatesens = 'no';
        bic = input(['Select bad ICs for',subj,':']);
        cfg.component = comp.label(bic);
        data_fix = ft_rejectcomponent(cfg, comp, f_data);
        while satis == 0
            close all
            comp = vy_ica(data_fix,lay);
            satis = input('Statisfied with ICA, yes = 1, not yet = 0:');
            if satis == 1, close('all'), break, end
            bic = input(['Select bad ICs for',subj,':']);
            cfg.component = comp.label(bic);
            data_fix = ft_rejectcomponent(cfg, comp, data_fix);
        end
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
        save(savepath, 'data_fix', '-v7.3');
    end
        
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
%     savepath = fullfile(outputdir,['tfr_ica',subj,'_',run]);
%     vy_tfr(data_fix,lay,savepath);
       
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
%     a_data = vy_ave(data_fix);
%     savepath = fullfile(outputdir,['gmean_',subj,'_',run]);
%     vy_ave_plot(a_data,lay,savepath)
    %     vy_ave_plot(a_data,lay,[])
   
    %% SPM source analysis
    data_clean = data_fix;
    
    switch subj
        case 'TQuality1'
            subj = 'C-117';
        case 'TQuality2'
            subj = 'C-120';
    end
    
    addpath(genpath(spm_path))
    cfg = [];
    cfg.toilim = [-0.4 1];
    eint_data = ft_redefinetrial(cfg, data_clean);
    
    mripath = fullfile(DestDirectory,'MRI','Nifti','T1',[subj,'_T1.nii']);
    name = [subj,'_',run];
    
    outputdir1 = fullfile(outputdir, 'spm_source');
    if exist(outputdir1, 'file') == 0
        mkdir(outputdir1);   % create a directory
    end
    
    %%
    if exist(mripath, 'file') == 2
        cd(outputdir1);
        vy_forward_spm_meg(datafile,eint_data,mripath,subj);
        cd(outputdir)
        
        %% from SPM (surface-based)
        disp('Source model ...');
        savepath_sourcemodel = ['sourcemodel1_',subj,'_',run,'.mat'];
%         if exist(savepath_sourcemodel, 'file') == 2
%             load(savepath_sourcemodel)
%         else
            addpath(genpath(spm_path))
            d = ['spm_source\m',subj];
            D = spm_eeg_load(d);
            datareg = D.inv{1}.datareg;
            mesh = spm_eeg_inv_transform_mesh(datareg.fromMNI*D.inv{1}.mesh.Affine, D.inv{1}.mesh);
            allmeshvert_mni = D.inv{1}.mesh.tess_mni.vert;
            
            hsfile = fullfile(datafolder1,'hs_file'); % headshape
            headshape = ft_read_headshape(hsfile);
            headshape = ft_convert_units(headshape, 'mm');
            
            % create the source model
            sourcemodel = [];
            sourcemodel.pos = mesh.tess_ctx.vert;
            sourcemodel.tri = mesh.tess_ctx.face;
            sourcemodel.inside = true(size(sourcemodel.pos,1),1);
            sourcemodel.unit = 'mm';
            
            % create the head model
            hdm                = [];
            hdm.bnd.pos        = mesh.tess_iskull.vert; % template's inner skull surface
            hdm.bnd.tri        = mesh.tess_iskull.face;
            hdm.type           = 'singleshell';
            hdm.unit           = 'mm';
            
            % compute the leadfield
            cfg             = [];
            cfg.senstype    = 'MEG';
            cfg.grid        = sourcemodel;
            cfg.headmodel   = hdm;
            cfg.channel     = {'MEG'};
            lf              = ft_prepare_leadfield(cfg, t_data.app);
            
            individual_grid = lf;
            individual_headmodel = hdm;
            
            %     figure;
            %     ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
            %     hold on;
            %     ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
            %     view ([-10 40 10])
            
            save(savepath_sourcemodel, 'individual_headmodel','individual_grid', 'allmeshvert_mni','-v7.3');
%         end
        %- anatomy inspections ( headmodel, mesh, alignment, ...)
        %  vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,sens, outputmridir)
        
        %%
        disp('wavelet analysis ...');
        savepath_wavelet = ['wavelet1_',subj,'_',run,'.mat'];
        if exist(savepath_wavelet, 'file') == 2
            load(savepath_wavelet)
        else
            foi = 1:1:30;
            w_data = vy_wavelet(data_clean, lay, savepath_wavelet,foi);
        end
        
        %%
        outputdir1 = fullfile(outputdir, 'all_sources');
        if exist(outputdir1, 'file') == 0, mkdir(outputdir1), end
        mtd = 'all_sources';
        d = ['.\spm_source\m',subj];
        D = spm_eeg_load(d);
        allmeshvert_mni = D.inv{1}.mesh.tess_mni.vert;
        cd(outputdir1);
        vy_source_all
%         vy_source_network
        %     vy_source_all2
        cd(outputdir)
        close all
        vy_init(cd_org)
    end
    
end



