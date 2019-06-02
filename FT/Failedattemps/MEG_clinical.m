clear; clc, close('all'); warning off

%%
restoredefaultpath

%% Initial settings, tools and data directory
addpath(genpath('./functions'));
addpath(genpath('./Data_file'));

cd_org = cd;
cd_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
cd (cd_tools)

%- Adding path
cfg = [];
cfg.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
allpath = vy_init(cfg);

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

%%
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
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

%-- excluding eyes-closed samples
% d = rdir([datadir,['/**/','sss','/*',tag,'*_closed/*_tsss.fif']]);
% clear datafolder_ex datafile_ex
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     datafolder_ex{i} = pathstr;
%     datafile_ex{i} = d(i).name;
% end
% datafolder_ex = datafolder_ex';
% disp(datafolder_ex)
%
% %%
% Acommon = intersect(datafile1,datafile_ex);
% datafile_sel = setxor(datafile1,Acommon);
% disp(datafile_sel)

%%
epoch_type = 'STI101';
% run = '1';

%%
disp('1: LCMV source')
disp('2: Network/Connectvity - Broadband');
disp('3: DICS Source, Beta');
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
% cfg.layout = 'neuromag306planar.lay';
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
    
    %%
    event = ft_read_event(datafile);
    clear val
    for i=1:length(event)
        val(i) = event(i).value;
    end
    val1 = unique(val);
    %     length(find(val == val1(1)))
    %     length(find(val == val1(2)))
    Evnt_IDs = min(val1);
    
    %%
%     cfg                         = [];
%     cfg.dataset                 = datafile;
% %     cfg.channel = {'MISC101'};
%     raw_data = ft_preprocessing(cfg);
%     
%     cfg = [];
%     cfg.viewmode = 'vertical';
%     cfg.continuous = 'no';
%     ft_databrowser(cfg,raw_data);
    
    %% Filteting, Event reading, Artifact rejecrtion
    savepath = fullfile(outputdir,[mtag,'_f_',subj,'.mat']);
    if exist(savepath, 'file') == 2
        load(savepath)
    else
        disp('preprocessing ...');
        switch method
            case 1
                f_data = vy_preprocess_beta (datafile,Evnt_IDs,epoch_type);
                f_data = vy_preprocess_gamma (datafile,Evnt_IDs,epoch_type);
            case {2,3}
                [f_data, ecg_data] = vy_preprocess(datafile,Evnt_IDs,epoch_type);
        end
        disp('preprocessing was completed');
        save(savepath, 'f_data', '-v7.3');
    end
    
    %% ICA cleaning
    savepath = fullfile(outputdir,[mtag,'_ica_',subj,'.mat']);
    saveflag = 1;
    ic_data = vy_ica_cleaning(f_data, [], lay, savepath, saveflag);
    
    %%
    disp('rej bad /channels/trials ...');
    savepath = fullfile(outputdir,[mtag,'_cl_',subj,'.mat']);
    saveflag = 1;
    cln_data = vy_rejectvisual(ic_data, lay, savepath, saveflag);
    
    %% data inspection
    %     cfg = [];
    %     cfg.viewmode = 'vertical';
    %     cfg.continuous = 'yes';
    %     ft_databrowser(cfg,cln_data);
    
    %% freq analysis (fft)
    savepath = fullfile(outputdir,[mtag,'_fft_',subj,'.mat']);
    saveflag = 1;
    vy_fft(cln_data, [2 40],1,savepath,saveflag);
    
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
            method = input('Eneter toi (e.g. [-0.3,0;1.5,2]): ');
%             toi = [-0.3,0;1.5,2];
            toi = [-0.3,0;0.5,0.8];           
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
    
    %% Data analysis
    switch analysis
        case 1
            %- Surface-based analysis
            vy_sourfacebasedsource
        case 2
            %- Volumetric-based analysis
            mridir = fullfile(datadir,subj,'brainstorm_db/anat');
            d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
            if ~isempty(d)
                sMRI1 = d.name;
                load(sMRI1); fid = SCS;
                mripfile = fullfile(mridir,'T1.nii');
                outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
                if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end               
                hsfile = datafile; % headshape
                [mri_realigned,individual_seg,individual_headmodel,headshape] = vy_mri_neuromag(mripfile,hsfile,fid,outputmridir,subj);                
                if exist('mri_realigned', 'var') == 1                   
                    disp('1: Low-res grid')
                    disp('2: High-res grid')
                    meshgrid = input('Eneter the mesh grid: ');
                    switch meshgrid
                        case 1
                            meshtag = 'lowres';
                        case 2
                            meshtag = 'highres';
                    end
                    
                    switch meshgrid
                        case 1
                            load temp_grid % low-res
                        case 2
                            load temp_grid_8mm % high-res
                            %                         temp = load('standard_sourcemodel3d8mm');
                    end
                    
                    %%
                    %                     cfg                        = [];
                    %                     cfg.grad                   = t_data.app.grad;
                    %                     cfg.headmodel              = individual_headmodel;
                    %                     cfg.reducerank             = 2;
                    %                     cfg.resolution             = 12;   % use a 3-D grid with a 1 cm resolution
                    %                     cfg.sourcemodel.unit       = 'mm';
                    %                     grid = ft_prepare_leadfield(cfg);                   
                    
                    %% Source modelling (warped with template)
                    [status] = ft_hastoolbox('SPM12', 1, 0);
                    cfg                 = [];
                    cfg.grid.warpmni    = 'yes';
                    cfg.spmversion     = 'SPM12';
                    cfg.grid.nonlinear  = 'yes';
                    cfg.grid.template   = template_grid;
                    cfg.mri             = mri_realigned;
                    cfg.grid.unit       = 'mm';
                    individual_grid     = ft_prepare_sourcemodel(cfg);
                    
                    %%
                    % make a figure with the template head model and dipole grid
                    %                 template_grid = ft_convert_units(template_grid, 'mm');
                    %                 figure
                    %                 hold on
                    %                 ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
                    %                 ft_plot_mesh(template_grid.pos(template_grid.inside,:));
                    
                    %% Anatomoy check!
                    anatomy_check_flag = 1;
                    saveflag = 2;
                    if anatomy_check_flag == 1
                        vy_mri_inspection(individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
                    end
                    close all
                    
                    %%
                    outputdir = fullfile(outputdir,mtag);
                    switch method
                        case 1
                            vy_source_lcmv_light
                            %                             vy_source_lcmv
                        case 2
                            %                             vy_network
                            vy_network_light1 % conn-network analysis
                        case 3
                            vy_source_dics
                        case 4
                            outputdir = fullfile(outdir,'ft_process',yttag, subj, tag);
                            outputdir1 = fullfile(outputdir, 'spm_source');
                            if exist(outputdir1, 'file') == 0
                                mkdir(outputdir1);   % create a directory
                            end
                            
                            addpath(genpath(spm_path))
                            cfg = [];
                            cfg.toilim = [-0.4 1];
                            eint_data = ft_redefinetrial(cfg, cln_data);
                            
                            if exist(mripfile, 'file') == 2
                                cd(outputdir1);
                                vy_forward_spm_meg(datafile,eint_data,mripfile,subj);
                                cd(outputdir)
                            end
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
end





