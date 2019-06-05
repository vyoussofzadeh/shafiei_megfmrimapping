clear; clc, close('all'); warning off

%%
restoredefaultpath

%% Initial settings
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));
cd_org = cd;
cd_tools = '/data/MEG/Vahab/Scripts/Vahab/Scripts/tools';

cd (cd_tools)
vy_init

datadir = '/data/MEG/Vahab/Scripts/Vahab/';

% datadir = '/home/vyoussof/Desktop/Vahab';

% outputdir1 = '/home/vyoussof/Desktop/Vahab/Processed_data'; % saving directory

cd(cd_org)

%%
clc;
% subj = 'FD';
subj = 'PB';
epoch_type = 'STI101';

run = '1';
mridir = fullfile(datadir,'raw_data',subj);

disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

%%
switch task
    case 1
        % % - Auditory definition naming
        tag = 'DFNAM';
        Evnt_IDs = 1; % questions
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
end

%%
d = rdir([mridir,['/*_',tag,'*'],'/*_tsss.fif']);
clear datafolder
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
end
datafolder = datafolder';
disp(datafolder)
datafile = d(i).name;

%%
% switch task
%     case 1
%         % % - Auditory definition naming
%         tag = 'DFNAM';
%         datafolder = fullfile(mridir,'Run09_DFNM_vertical');
%         datafile = fullfile(datafolder,'Run09_DFNM_vertical_raw_tsss.fif');
%         Evnt_IDs = 1; % questions
%     case 2
%         % - Visual picture naming
%         tag = 'PicNam';
%         datafolder = fullfile(mridir,'Run08_PN_vertical');
%         datafile = fullfile(datafolder,'Run08_PN_vertical_raw_tsss.fif');
%         Evnt_IDs = [2,3]; % 3: images, 2: scrambled images
% end

%%
disp('1: Beta source')
disp('2: Network-connectvity analysis');
disp('3: Dics bf - beta');
method = input('Eneter the method: ');
switch method
    case 1
        mtag = 'beta';
    case 2
        mtag = 'conn';
    case 3
        mtag = 'dics';
end

%%
disp('1: Low-res grid')
disp('2: High-res grid')
meshgrid = input('Eneter the mesh grid: ');
switch meshgrid
    case 1
        meshtag = 'lowres';
    case 2
        meshtag = 'highres';
end

%%
outputdir = fullfile(datadir,'processed_data',subj, tag);
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end
cd(outputdir)

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%%
for i = 1:size(datafolder,1)
    
    disp([datafolder,' is analysing']),
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outputdir,[mtag,'_f_',subj,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        switch method
            case 1
                f_data = vy_preprocess_beta (datafile,Evnt_IDs,epoch_type);
            case {2,3}
                f_data = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        end
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
    %% ICA cleaning
    savepath = fullfile(outputdir,[mtag,'_ica_',subj,'.mat']);
    saveflag = 1;
    ic_data = vy_ica_cleaning(f_data, lay, savepath, saveflag);
    
    %%
    disp('rej bad /channels/trials ...');
    savepath = fullfile(outputdir,[mtag,'_cl_',subj,'.mat']);
    saveflag = 1;
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
    savepath = fullfile(outputdir,[mtag,'_fft_',subj,'.mat']);
    saveflag = 1;
    vy_fft(cln_data, [2 40],1,savepath,saveflag);
    
    % pausue
    %% freq analysis (tfr)
    savepath = fullfile(outputdir,[mtag,'_tft_',subj,'.mat']);
    vy_tfr(cln_data,lay,savepath);
    %     vy_tfr(data_clean,lay,[]);
    
    %% elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    % data.sens = sens;
    
    %% Epoching
    switch task
        case 1 %'DefNam'
            toi = [-0.3,0;1.5,2];
        case 2 %'PicNam'
            %             toi = [-0.3,0;0.7,1];
            toi = [-0.2,0;0.5,0.8];
    end
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
    savepath = fullfile(outputdir,[mtag,'_gmean_',subj,'.mat']);
    vy_ave_plot(a_data,lay,savepath)
    
    %% mri anatomy
    mripath = fullfile(mridir,'T1.nii');
    outputmridir = fullfile(datadir,'processed_data',subj,'anat'); % output dir
    
    if exist(outputmridir, 'file') == 0
        mkdir(outputmridir);   %create the directory
    end
    hsfile = datafile; % headshape
    cd '/data/MEG/Vahab/Scripts/Vahab/Test_MEG2/bednar_peggy/brainstorm_db/anat/bednar_peggy'
    load subjectimage_T1.mat
    fid = SCS;
    
    if exist(mripath, 'file') == 2
        
        [mri_realigned,individual_seg, individual_headmodel,headshape] = vy_mri_neuromag(mripath,hsfile,fid,outputmridir,subj);
        if exist('mri_realigned', 'var') == 1
            
            switch meshgrid
                case 1
                    load temp_grid % low-res
                case 2
                    load temp_grid_8mm % high-res
            end
            
            %% Source modelling (warped with template)
            cfg                 = [];
            cfg.grid.warpmni    = 'yes';
            cfg.grid.nonlinear  = 'yes';
            cfg.grid.template   = template_grid;
            cfg.mri             = mri_realigned;
            cfg.grid.unit       = 'mm';
            individual_grid     = ft_prepare_sourcemodel(cfg);
            
            %% Anatomoy check!
            anatomy_check_flag = 1;
            saveflag = 1;
            if anatomy_check_flag == 1
                vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
            end
            close all
            
            %%
            outputdir = fullfile(outputdir,mtag);
            switch method
                case 1
                    %                     vy_source_lcmv_light
                    vy_source_lcmv
                case 2
                    %             vy_network
                    vy_network_light1 % conn-network analysis
                case 3
                    vy_source_dics
            end
            
            %%
            clc
            close all
            disp([datafile,' ,was completed'])
            
        else
            disp([datafile,' ,was not completed'])
        end
    end
end




