function varargout = vy_ft_sourceanalysis_dics_vol( varargin )
% PROCESS_FT_SOURCEANALYSIS Call FieldTrip function ft_sourceanalysis

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2016-2017

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

% ===== PROCESS =====
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 325;
sProcess.Description = 'http://www.fieldtriptoolbox.org/tutorial/minimumnormestimate';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Label: Warning
sProcess.options.label1.Comment = '<B>Warning</B>: Process under development.<BR><BR>';
sProcess.options.label1.Type    = 'label';
% Option: Inverse method
sProcess.options.method.Comment = 'Inverse method:';
sProcess.options.method.Type    = 'combobox_label';
sProcess.options.method.Value   = {'mne', {'LCMV beamformer', 'SAM beamformer', 'DICS beamformer', 'MNE', 'sLORETA', 'eLORETA', 'MUSIC', 'PCC', 'Residual variance'; ...
    'lcmv',            'sam',            'dics',            'mne', 'sloreta', 'eloreta', 'music', 'pcc', 'rv'}};
% Option: Sensors selection
sProcess.options.sensortype.Comment = 'Sensor type:';
sProcess.options.sensortype.Type    = 'combobox_label';
sProcess.options.sensortype.Value   = {'MEG', {'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'; ...
    'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'}};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};
% Initialize fieldtrip
bst_ft_init();

% ===== GET OPTIONS =====
% Inverse options
Method   = sProcess.options.method.Value;
Modality = sProcess.options.sensortype.Value{1};
% Get unique channel files
AllChannelFiles = unique({sInputs.ChannelFile});
% Progress bar
bst_progress('start', 'ft_sourceanalysis', 'Loading input files...', 0, 2*length(sInputs));

% ===== LOOP ON FOLDERS =====
for iChanFile = 1:length(AllChannelFiles)
    bst_progress('text', 'Loading input files...');
    % Get the study
    [sStudyChan, iStudyChan] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
    % Error if there is no head model available
    if isempty(sStudyChan.iHeadModel)
        bst_report('Error', sProcess, [], ['No head model available in folder: ' bst_fileparts(sStudyChan.FileName)]);
        continue;
    elseif isempty(sStudyChan.NoiseCov) || isempty(sStudyChan.NoiseCov(1).FileName)
        bst_report('Error', sProcess, [], ['No noise covariance matrix available in folder: ' bst_fileparts(sStudyChan.FileName)]);
        continue;
    end
    % Load channel file
    ChannelMat = in_bst_channel(AllChannelFiles{iChanFile});
    % Get selected sensors
    iChannels = channel_find(ChannelMat.Channel, Modality);
    if isempty(iChannels)
        bst_report('Error', sProcess, sInput, ['Channels "' Modality '" not found in channel file.']);
        return;
    end
    % Load head model
    HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
    HeadModelMat = in_bst_headmodel(HeadModelFile);
    % Load data covariance matrix
    NoiseCovFile = sStudyChan.NoiseCov(1).FileName;
    NoiseCovMat = load(file_fullpath(NoiseCovFile));
    %%% DATA OR NOISE COVARIANCE ????
    
    %%
    %     set(0,'DefaultFigureWindowStyle','docked')
    %
    %     cd('/MEG_data/Vahab/Github/MCW-MEGlab/FT');
    %     restoredefaultpath
    %     cd_org = cd;
    %     addpath(genpath(cd_org));
    %
    %     %- Input dir
    %     indir = '/MEG_data/epilepsy';
    %     %- Output dir
    %     outdir = '/MEG_data/Vahab/Processed_data';
    %
    %     %- Adding path
    %     cfg_init = [];
    %     cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
    %     [allpath, atlas] = vy_init(cfg_init);
    
    %%
    Index = strfind(HeadModelMat.SurfaceFile, '/');
    subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
    indir = '/MEG_data/epilepsy';
    bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
    bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
    bsdir = fullfile(indir,subj,'brainstorm_db');
    cd(bsdir)
    
    %%
    %     close all
    %     %- Adding path
    %     cfg_init = [];
    %     cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
    %     [allpath, atlas] = vy_init(cfg_init);
    %%
    %     brainstorm stop
    % rmpath(genpath('/usr/local/MATLAB_Tools/brainstorm3'));
    
    %% Reading trials
    % ===== LOOP ON DATA FILES =====
    % Get data files for this channel file
    iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
    % Loop on data files
    for iInput = 1:length(iChanInputs)
        
        % === LOAD DATA ===
        % Load data
        DataFile = sInputs(iChanInputs(iInput)).FileName;
        DataMat = in_bst_data(DataFile);
        iStudyData = sInputs(iChanInputs(iInput)).iStudy;
        % Remove bad channels
        iBadChan = find(DataMat.ChannelFlag == -1);
        iChannelsData = setdiff(iChannels, iBadChan);
        
        if isempty(iChannelsData)
            bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
            return;
        end
        
        trl{iInput} = DataMat.F;
        timee {iInput} = DataMat.Time;
        
    end
    
    ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
    ftData1 = [];
    ftData1.label = ftData.label;
    ftData1.grad = ftData.grad;
    ftData1.dimord = ftData.dimord;
    ftData1.trial = trl;
    ftData1.time = timee;
    
    %%
    Index = strfind(DataFile, '/');
    saveid = DataFile(Index(end)+1:end-4);
    savepath = fullfile(['dics_',saveid]);
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    
    %% Selecting freq of interest
    if exist(fullfile(savepath,'tfr.mat'),'file')~= 2
        % do tfr-decomposition
        cfg = [];
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.method     = 'mtmconvol';
        cfg.taper      = 'hanning';
        % cfg.taper      = 'dpss';
        cfg.foi        = 1:3:40;
        cfg.keeptrials = 'yes';
        cfg.t_ftimwin  = 3./cfg.foi;
        % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.tapsmofrq  = 0.8 *cfg.foi;
        cfg.toi        = -0.5:0.05:2;
        tfr        = ft_freqanalysis(cfg, ftData1);
        save(fullfile(savepath,'tfr.mat'),'tfr');
    else
        load(fullfile(savepath,'tfr.mat'));
    end
    
    cfg = [];
    cfg.savepath = [];
    cfg.savefile = [];
    [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
    disp(['time_of_interest:',num2str(time_of_interest),'sec']);
    disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
    L = 0.3;
    
    %% Selecting time of interest
    datain = ftData1;
    %     %--setting baseline interval
    %     toi(1,:) = [-0.3,0];
    %
    %     %-- setting the post-stim interval
    %     disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
    %     %     disp('Yes: 1, No: 2');
    %     %     tfa = input('Is it OK to proceed?');
    %     %     disp('the following time was selected');
    %     %        if (time_of_interest-L) > 0.4 && (time_of_interest+L) < 2.3
    %     if (time_of_interest-L) >= 0.1 && (time_of_interest+L) < 2.3
    %         toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
    %     else
    %         switch task
    %             case 1
    %                 %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
    %                 toi(2,:) = [0.8,1.5]; % Best of DN
    %             case 2
    %                 toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
    %                 toi(2,:) = [0.4,1.6]; % Best for PN, left IFG
    %                 toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
    %                 toi(2,:) = [1,1.6]; % Best for PN, left IFG
    %             case 3
    %                 toi(1,:) = [-0.25,0.25];
    %                 toi(2,:) = [0.75,1.25]; % Best for PN, left IFG
    %         end
    %     end
    %     disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
    %     disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
    %     % Run_time
    %     ep_data = vy_epoch(datain, toi);
    %
    %     %- Appending data
    %     cfg = [];
    %     ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %--setting baseline interval
    toi(1,:) = [-0.3,0];
    
    %-- setting the post-stim interval
    disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
    %     disp('Yes: 1, No: 2');
    %     tfa = input('Is it OK to proceed?');
    %     disp('the following time was selected');
    %        if (time_of_interest-L) > 0.4 && (time_of_interest+L) < 2.3
    if (time_of_interest-L) >= 0.1 && (time_of_interest+L) < 2.3
        toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
    else
        switch task
            case 1
                %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
                toi(2,:) = [0.8,1.5]; % Best of DN
            case 2
                toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
                toi(2,:) = [0.4,1.6]; % Best for PN, left IFG
                toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
                toi(2,:) = [1,1.6]; % Best for PN, left IFG
            case 3
                toi(1,:) = [-0.25,0.25];
                toi(2,:) = [0.75,1.25]; % Best for PN, left IFG
        end
    end
    disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
    disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
    Run_time
    
    
    %%
%     cfg = [];
%     cfg.layout = 'neuromag306mag.lay';
%     lay = ft_prepare_layout(cfg);
%     Run_grandmean
%     
%     %%
%     %     cfg_main = [];
%     %     cfg_main.outputdir = [];
%     %     cfg_main.freq_of_interest = freq_of_interest;
%     %     cfg_main.sens = datain.grad;
%     %     Run_fft_4dics
%     
%     %%
%     
%     
%     %% Volumetric-based analysis
%     Index = strfind(HeadModelMat.SurfaceFile, '/');
%     subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
%     
%     outdir = '/MEG_data/Vahab/Processed_data';
%     outputmridir = fullfile(outdir,'ft_process', subj,'anat'); % output dir
%     d = rdir(fullfile(bsanatdir,subj,'subjectimage*.mat'));
%     clear fid
%     if ~isempty(d)
%         sMRI1 = d.name;
%         load(sMRI1);
%         fid.SCS = SCS;
%         fid.NCS = NCS;
%         mripfile = fullfile(bsanatdir,'T1.nii');
%         cfg = [];
%         cfg.megdata = ep_data.app;
%         cfg.mripfile = mripfile;
%         cfg.hsfile = []; % headshape;
%         cfg.fid = fid;
%         cfg.outputmridir = outputmridir;
%         cfg.subj = subj;
%         cfg.plotflag = 2;
%         outanat = vy_mri_neuromag2(cfg);
%     end
%     
%     %% Choosing mesh
%     if flag.warping == 1
%         switch meshgrid
%             case 1
%                 meshtag = 'lowres';
%                 %         load('standard_sourcemodel3d10mm');
%                 load temp_grid % low-res
%                 template_grid = ft_convert_units(template_grid, 'mm');
%                 individual_grid = outanat.individual_grid_10mm;
%             case 2
%                 meshtag = 'highres';
%                 %         load('standard_sourcemodel3d8mm');
%                 load temp_grid_8mm % high-res
%                 individual_grid = outanat.individual_grid_8mm;
%         end
%     else
%         switch meshgrid
%             case 1
%                 individual_grid = outanat.individual_grid_10mm;
%             case 2
%                 individual_grid = outanat.individual_grid_8mm;
%         end
%     end
%     
%     %% Anatomoy check!
%     saveflag = 2; flag.anatomycheck = 1;
%     if flag.anatomycheck == 1
%         cfg = [];
%         cfg.saveflag = saveflag;
%         cfg.headmodel = outanat.individual_headmodel;
%         cfg.leadfield = individual_grid;
%         cfg.mri_realigned  = outanat.mri_realigned;
%         cfg.headshape = outanat.headshape;
%         cfg.outputmridir = outputmridir;
%         cfg.mtd = 'vol';
%         vy_mri_inspection(cfg);
%         %     vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
%     end
%     
%     %%
%     ft_path =  '/MEG_data/Vahab/Github/tools/ft_packages/fieldtrip_20190419';
%     template_mri = ft_read_mri(fullfile(ft_path,'template/anatomy','single_subj_T1.nii')); %
    
    %%
    %     cfg = [];
    %     cfg.grid = individual_grid;
    %     cfg.allpath = allpath;
    %     cfg.freq_of_interest  = freq_of_interest; % Hz
    %     cfg.headmodel = outanat.individual_headmodel;
    %     cfg.sens = sens;
    %     cfg.mtag = mtag;
    %     cfg.subj = subj;
    %     cfg.toi = toi;
    %     cfg.outputdir = outd.vol;
    %     if flag.warping ==1
    %         cfg.template_grid = template_grid;
    %     end
    %     cfg.template_mri = template_mri;
    %     cfg.savedata = fullfile(outd.vol,[mtag,'_',subj]);
    %     cfg.flag = flag;
    %     vy_source_dics(cfg, ep_data);
    
    %     %% Leadfield estimation
    %     cfg                 = [];
    %     cfg.grad            = ep_data.all.grad;
    %     cfg.headmodel       = outanat.individual_headmodel;
    %     cfg.reducerank      = 2;
    %     cfg.channel         = {'MEG'};
    %     cfg.resolution = 10;   % use a 3-D grid with a 1 cm resolution
    %     cfg.sourcemodel.unit       = 'mm';
    %     individual_grid_10mm = ft_prepare_leadfield(cfg);
    %
    %     figure;
    %     ft_plot_vol(outanat.individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    %     hold on;
    %     ft_plot_headshape(outanat.headshape);
    %     ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
    %     view ([0 90])
    %
    %     %%
    %     mri_idiv = ft_read_mri(fullfile(bsanatdir,'T1.nii')); %
    
    
    %%
%     load('/MEG_data/Vahab/Processed_data/ft_process/bednar_peggy/DFN/ep_data.mat')
    
    %%
%     individual_grid1 = ft_transform_geometry(ChannelMat.TransfMeg{2},individual_grid);
%     template_grid1 = ft_transform_geometry(ChannelMat.TransfMeg{2},template_grid);
%     individual_headmodel1 = ft_transform_geometry(ChannelMat.TransfMeg{2},outanat.individual_headmodel);
%     template_mri1 = ft_transform_geometry(ChannelMat.TransfMeg{2},template_mri);
%     
%     sourcemodel_mm = ft_convert_units(sourcemodel,'mm');
%     
%     
%     figure;
%     % ft_plot_vol(individual_headmodel_surf, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     ft_plot_vol(individual_headmodel1 , 'facecolor', 'cortex', 'edgecolor', 'none');     alpha 0.5;
%     hold on;
%     ft_plot_headshape(sourcemodel_mm);
%     ft_plot_mesh(individual_grid1.pos(individual_grid1.inside, :));
%     view ([-10 40 10])
    
%     cfg = [];
%     cfg.grid = individual_grid1;
%     %     cfg.allpath = allpath;
%     cfg.freq_of_interest  = freq_of_interest; % Hz
%     cfg.headmodel = individual_headmodel1;
%     cfg.sens = datain.grad;
%     cfg.mtag = 'dics';
%     cfg.subj = subj;
%     cfg.toi = toi;
%     cfg.outputdir = savepath;
%     if flag.warping ==1
%         cfg.template_grid = template_grid1;
%     end
%     cfg.template_mri = template_mri1;
%     cfg.savedata = fullfile(['dics_',subj]);
%     cfg.flag = flag;
%     vy_source_dics(cfg, ep_data);

cfg_main = [];
cfg_main.sens = datain.grad;
cfg_main.outputdir = savepath;
cfg_main.freq_of_interest  = freq_of_interest; % Hz
Run_fft_4dics
%%
% cfg = [];
% cfg.headmodel = individual_headmodel1;
% % cfg.sourcemodel = cfg_main.sourcemodel;
% cfg.grid = individual_grid1;
% cfg.mtag = 'dics';
% s_data_dics = vy_source_freq(cfg, f_data);
% % s_data_dics = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, cfg_main.mtag);
% 
% %%
% cfg = [];
% cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x1+x2)';
% source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
% source_diff_dics.pos     = template_grid.pos;
% source_diff_dics.dim     = template_grid.dim;
% source_diff_dics.inside  = template_grid.inside;
% source_diff_dics.pow(source_diff_dics.pow>0)=0;
% 
% %%
% cfg = [];
% cfg.mask = 'pow';
% cfg.loc = 'min';
% cfg.template = template_mri;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     source_int_dics = vy_source_plot(cfg, source_diff_dics);
    
%     %%
%     cfg = [];
%     cfg.parameter = 'pow';
%     % cfg.interpmethod = 'sphere_avg';
%     cfg.interpmethod = 'smudge';
%     cfg.coordsys     = 'mni';
%     test = ft_sourceinterpolate(cfg, source_diff_dics, template_mri);
%     
%     
%     
%     mesh1 = [];
%     mesh1.unit = 'm';
%     mesh1.tri = sourcemodel.tri;
%     mesh1.sulc = sourcemodel.sulc;
%     mesh1.curve = sourcemodel.curv;
%     mesh1.pos = sourcemodel.pos;
%     
%     a.mesh
%     
%     cfg                = [];
%     cfg.method         = 'ortho';
%     cfg.funparameter   = 'pow';
%     cfg.funcolorlim    = 'maxabs';
%     cfg.opacitymap     = 'rampup';
%     cfg.crosshair      = 'no';
%     cfg.camlight       = 'no';
%     cfg.funcolormap    =  brewermap(256, '*RdYlBu');
%     cfg.projthresh     = 0.6;   
%     cfg.method = 'surface';
%     cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% %     cfg.surfinflated   = sourcemodel_mm;
% %     cfg.surfinflated   = 'surface_inflated_right.mat';
%     ft_sourceplot(cfg, test);
%     view([-90 0]); camlight; material dull;
%     
%     %%
%     maskval = ones(size(tmpdata.pow,1),1);
%     figure,
%     ft_plot_mesh(sourcemodel_mm, 'edgecolor', 'none', 'facecolor', [], 'vertexcolor', 'curv');
%     ft_plot_mesh(sourcemodel_mm, 'edgecolor', 'none', 'facecolor', [], ...
%         'vertexcolor', tmpdata.pow, 'facealpha', maskval, 'clim', [-0.5372    0.5372], ...
%         'alphalim', [ 0, 1], 'alphamap', ft_getopt(cfg, 'opacitymap',    'auto'), 'colormap', ...
%         [], 'maskstyle', 'opacity');    
%     
%     %%
%     tmpcfg = [];
%     tmpcfg.parameter = 'pow';
%     tmpcfg.interpmethod = 'sphere_avg';
% %     tmpcfg.distmat      = [];
%     tmpcfg.sphereradius = 3;
% %     tmpcfg.projvec      = 1;
%     tmpcfg.projcomb     = 'mean';
% %     tmpcfg.projweight   = 1;
% %     tmpcfg.projthresh   = 0;
%     tmpdata             = ft_sourceinterpolate(tmpcfg, source_diff_dics, sourcemodel_mm);
%     
%      
%     figure
%     m = tmpdata.pow;
%     bnd.pnt = tmpdata.pos;
%     bnd.tri = tmpdata.tri;
%     ft_plot_mesh(bnd, 'vertexcolor', abs(m));
%     colorbar
%     
%     %%
%     surf = sourcemodel;
%     ft_plot_mesh(surf, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;

    
    %%
    %     figure;
    %     ft_plot_vol(outanat.individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    %     hold on;
    %     %     ft_plot_headshape(cfg.headshape);
    %     ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
    %     view ([0 90])
    %
    %     %     a = [];
    %     %     a.pos = ChannelMat.HeadPoints.Loc';
    %     %     a.label = ChannelMat.HeadPoints.Label;
    %     %     a = ft_convert_units(a,'mm');
    %     %     ft_plot_headshape(a);
    %
    %     ChannelMat1 = ChannelMat;
    %     ChannelMat1.unit = 'm';
    %     ChannelMat1 = ft_convert_units(ChannelMat1,'mm');
    %     HeadPoints = ft_transform_geometry(ChannelMat.TransfMeg{2}/outanat.mri_realigned.transform,ChannelMat1.HeadPoints);
    %     a = [];
    %     a.pos = HeadPoints.Loc';
    %     a.label = HeadPoints.Label;
    %     a = ft_convert_units(a,'mm');
    %     ft_plot_headshape(a);
    %
    %
    %     %%
    %     ep_data1 = ep_data;
    %     ep_data1.all.grad = ft_transform_geometry(ChannelMat.TransfMeg{1},ep_data1.all.grad);
    %     ep_data1.app.grad = ft_transform_geometry(ChannelMat.TransfMeg{1},ep_data1.app.grad);
    %     ep_data1.bsl.grad = ft_transform_geometry(ChannelMat.TransfMeg{1},ep_data1.bsl.grad);
    %     ep_data1.pst.grad = ft_transform_geometry(ChannelMat.TransfMeg{1},ep_data1.pst.grad);
    %
    %     ep_data1.all.grad = ft_transform_geometry(outanat.mri_realigned.transform,ep_data1.all.grad);
    %     ep_data1.app.grad = ft_transform_geometry(outanat.mri_realigned.transform,ep_data1.app.grad);
    %     ep_data1.bsl.grad = ft_transform_geometry(outanat.mri_realigned.transform,ep_data1.bsl.grad);
    %     ep_data1.pst.grad = ft_transform_geometry(outanat.mri_realigned.transform,ep_data1.pst.grad);
    %
    %     cfg = [];
    %     cfg.grid = individual_grid;
    %     %     cfg.allpath = allpath;
    %     cfg.freq_of_interest  = freq_of_interest; % Hz
    %     cfg.headmodel = outanat.individual_headmodel;
    %     cfg.sens = datain.grad;
    %     cfg.mtag = 'dics';
    %     cfg.subj = subj;
    %     cfg.toi = toi;
    %     cfg.outputdir = savepath;
    %     if flag.warping ==1
    %         cfg.template_grid = template_grid;
    %     end
    %     cfg.template_mri = template_mri1;
    %     cfg.savedata = fullfile(['dics_',subj]);
    %     cfg.flag = flag;
    %     vy_source_dics(cfg, ep_data1);
    
    
    %% Now on the surface (DICS on the surface)
%     sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
%     
%     % - Leadfield
%     if exist(fullfile(savepath,'leadfield.mat'),'file')~= 2
%         [individual_headmodel_surf, individual_grid_surf] = vy_bs2ft_headmodel(ep_data.all, sourcemodel);
%         save(fullfile(savepath,'leadfield.mat'),'individual_headmodel_surf', 'individual_grid_surf');
%     else
%         load(fullfile(savepath,'leadfield.mat'));
%     end
%     
%     sourcemodel1 = ft_convert_units(sourcemodel, 'mm');
%     individual_headmodel2 = ft_convert_units(individual_headmodel1, 'm');
%     
%     
%     figure;
%     % ft_plot_vol(individual_headmodel_surf, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     ft_plot_vol(individual_headmodel2 , 'facecolor', 'cortex', 'edgecolor', 'none');     alpha 0.5;
%     hold on;
%     ft_plot_headshape(sourcemodel);
%     ft_plot_mesh(individual_grid_surf.pos(individual_grid_surf.inside, :));
%     view ([-10 40 10])
       
    %%
       sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
        [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
     
        
             cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '0.1%';
        cfg.sourcemodel  = ftLeadfield;
        cfg.frequency    = f_data.app.freq;
        cfg.headmodel = ftHeadmodel;
        cfg.dics.keepfilter = 'yes';
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        sourceavg = ft_sourceanalysis(cfg, f_data.app);
        
        cfg = [];
        cfg.method = 'dics';
        cfg.dics.lambda = '0.1%';
        cfg.sourcemodel        = ftLeadfield;
        cfg.sourcemodel.filter = sourceavg.avg.filter;
        cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.headmodel = ftHeadmodel;
        s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
        s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
        
        %
        cfg = [];
        cfg.parameter = 'pow';
        %     cfg.operation = '(x1-x2)/(x1+x2)';
        cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
        source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
        source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
        source_diff_dics.pow(source_diff_dics.pow>0)=0;
        source_diff_dics.pow = (source_diff_dics.pow);
        
        figure
        m = source_diff_dics.pow;
        bnd.pnt = sourcemodel.pos;
        bnd.tri = sourcemodel.tri;
        ft_plot_mesh(bnd, 'vertexcolor', abs(m));
        colorbar
      
    %     cfg = [];
    %     cfg.headmodel = individual_headmodel2;
    %     % cfg.sourcemodel = cfg_main.sourcemodel;
    %     cfg.grid = individual_grid_surf;
    %     cfg.mtag = 'dics';
    %     s_data_dics = vy_source_freq(cfg, f_data);
    %
    %     %%
    %     cfg = [];
    %     cfg.parameter = 'pow';
    %     cfg.operation = '(x1-x2)/(x1+x2)';
    %     source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
    %     source_diff_dics.pos     = template_grid.pos;
    %     source_diff_dics.dim     = template_grid.dim;
    %     source_diff_dics.inside  = template_grid.inside;
    %     source_diff_dics.pow(source_diff_dics.pow>0)=0;
    
    %% DICS surface
%     cfg = [];
%     cfg.method = 'dics';
%     cfg.dics.lambda = '1%';
%     cfg.sourcemodel  = individual_grid_surf;
%     cfg.frequency    = f_data.app.freq;
%     cfg.headmodel = individual_headmodel2;
%     cfg.dics.keepfilter = 'yes';
%     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
%     sourceavg = ft_sourceanalysis(cfg, f_data.app);
%     
%     cfg = [];
%     cfg.method = 'dics';
%     cfg.dics.lambda = '1%';
%     cfg.sourcemodel        = individual_grid_surf;
%     cfg.sourcemodel.filter = sourceavg.avg.filter;
%     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
%     cfg.headmodel = individual_headmodel2;
%     s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
%     s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
%     
%     %
%     cfg = [];
%     cfg.parameter = 'pow';
%     cfg.operation = '(x1-x2)/(x1+x2)';
%     %     cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
%     source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
%     source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
%     source_diff_dics.pow(source_diff_dics.pow>0)=0;
%     source_diff_dics.pow = (source_diff_dics.pow);
%     
%     %%
%     figure
%     m = source_diff_dics.pow;
%     bnd.pnt = individual_grid_surf.pos;
%     bnd.tri = individual_grid_surf.tri;
%     ft_plot_mesh(bnd, 'vertexcolor', abs(m));
%     colorbar
%     %     hold on
    %     ft_plot_vol(individual_headmodel2 , 'facecolor', 'cortex', 'edgecolor', 'none');     alpha 0.5;
    
    
    %% Surface-based source analysis, from BS
%     [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
%     
%     figure;
%     % ft_plot_vol(individual_headmodel_surf, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     ft_plot_vol(ftHeadmodel , 'facecolor', 'cortex', 'edgecolor', 'none');     alpha 0.5;
%     hold on;
%     ft_plot_headshape(sourcemodel);
%     ft_plot_mesh(ftLeadfield.pos(ftLeadfield.inside, :));
%     view ([-10 40 10])
%     
%     cfg = [];
%     cfg.method = 'dics';
%     cfg.dics.lambda = '0.1%';
%     cfg.sourcemodel  = ftLeadfield;
%     cfg.frequency    = f_data.app.freq;
%     cfg.headmodel = ftHeadmodel;
%     cfg.dics.keepfilter = 'yes';
%     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
%     sourceavg = ft_sourceanalysis(cfg, f_data.app);
%     
%     cfg = [];
%     cfg.method = 'dics';
%     cfg.dics.lambda = '0.1%';
%     cfg.sourcemodel        = ftLeadfield;
%     cfg.sourcemodel.filter = sourceavg.avg.filter;
%     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
%     cfg.headmodel = ftHeadmodel;
%     s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
%     s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
%     
%     %
%     cfg = [];
%     cfg.parameter = 'pow';
% %     cfg.operation = '(x1-x2)/(x1+x2)';
%         cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
%     source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
%     source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
%     source_diff_dics.pow(source_diff_dics.pow>0)=0;
%     source_diff_dics.pow = (source_diff_dics.pow);
%     
%     
%     figure
%     m = source_diff_dics.pow;
%     bnd.pnt = individual_grid_surf.pos;
%     bnd.tri = individual_grid_surf.tri;
%     ft_plot_mesh(bnd, 'vertexcolor', abs(m));
%     colorbar
%     
%     %%
%     cfg = [];
%     cfg.mask = 'pow';
%     cfg.loc = 'min';
%     cfg.template = template_mri;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     source_int_dics = vy_source_plot(cfg, source_diff_dics);
%     set(gcf,'name',cfg_main.subj,'numbertitle','off')
%     
%     clear savepath
%     savepath{1} = fullfile(outputdir_dics,[num2str(f),'Hz','_2_',cfg_main.subj]);
%     savepath{2} = fullfile(outputdir_dics,[num2str(f),'Hz','_3_',cfg_main.subj]);
%     
%     cfg = [];
%     cfg.subj = subj;
%     cfg.mask = 'pow';
%     cfg.thre = 0.6;
%     cfg.savepath = savepath;
%     vy_mapvisualisation(cfg, source_int_dics);
%     % vy_mapvisualisation(source_int_dics,cfg.mask,0.6, []);
    
    %%
%     cfg            = [];
%     cfg.downsample = 2;
%     cfg.parameter = 'pow';
%     source_diff_dics_int  = ft_sourceinterpolate(cfg, source_diff_dics , outanat.mri_realigned);
%     
%     %%
%     cfg              = [];
%     cfg.method       = 'ortho';
%     cfg.funparameter = 'pow';
%     figure
%     ft_sourceplot(cfg,source_diff_dics_int);
%     
%     %%
%     cfg = [];
%     source_int_norm = ft_volumenormalise(cfg, source_diff_dics_int);
%     
%     %%
%     cfg = [];
%     cfg.subj = subj;
%     cfg.mask = 'pow';
%     cfg.thre = 0;
%     cfg.savepath = [];
%     vy_mapvisualisation(cfg, source_int_norm);
%     
%     %%
%     source_int1 = source_diff_dics_int;
%     tmp = abs(source_int1.pow);
%     tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
%     source_int1.pow = tmp;
%     
%     cfg = [];
%     % cfg.method        = 'ortho';
%     cfg.method        = 'slice';
%     cfg.funparameter  = 'pow';
%     %     cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = [0.1 1];
%     cfg.opacitylim    = [0.1 1];
%     cfg.nslices       = 16;
%     cfg.slicerange    = [10,60];
%     cfg.slicedim      = [3];
%     cfg.opacitymap    = 'rampup';
%     cfg.funcolormap =  brewermap(256, '*RdYlBu');
%     ft_sourceplot(cfg, source_int1);
%     
%     %%
%     cfg = [];
%     cfg.mask = 'pow';
%     cfg.loc = 'min';
%     cfg.template = mri_idiv;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     source_int_dics = vy_source_plot(cfg, source_int1);
%     
    %%
    %     sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
    %
    %     %% Leadfield
    %     if exist(fullfile(savepath,'leadfield.mat'),'file')~= 2
    %         [individual_headmodel, individual_grid] = vy_bs2ft_headmodel(ep_data.all, sourcemodel);
    %         save(fullfile(savepath,'leadfield.mat'),'individual_headmodel', 'individual_grid');
    %     else
    %         load(fullfile(savepath,'leadfield.mat'));
    %     end
    %
    %     %% DICS surface
    %     data = f_data;
    %     cfg = [];
    %     cfg.method = 'dics';
    %     cfg.dics.lambda = '1%';
    %     cfg.sourcemodel  = individual_grid;
    %     cfg.frequency    = data.app.freq;
    %     cfg.headmodel = individual_headmodel;
    %     cfg.dics.keepfilter = 'yes';
    %     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    %     sourceavg = ft_sourceanalysis(cfg, data.app);
    %
    %     cfg = [];
    %     cfg.method = 'dics';
    %     cfg.dics.lambda = '1%';
    %     cfg.sourcemodel        = individual_grid;
    %     cfg.sourcemodel.filter = sourceavg.avg.filter;
    %     cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    %     cfg.headmodel = individual_headmodel;
    %     s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
    %     s_data.pst      = ft_sourceanalysis(cfg, data.pst);
    %
    %     %
    %     cfg = [];
    %     cfg.parameter = 'pow';
    %     cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
    %     source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
    %     source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
    %     source_diff_dics.pow(source_diff_dics.pow>0)=0;
    %     source_diff_dics.pow = abs(source_diff_dics.pow);
    
    %% DICS volume
    
    
    %%
%     figure
%     m = source_diff_dics.pow;
%     bnd.pnt = individual_grid.pos;
%     bnd.tri = individual_grid.tri;
%     ft_plot_mesh(bnd, 'vertexcolor', sqrt(m));
%     colorbar
    
    %%
    
    % === CREATE OUTPUT STRUCTURE ===
    bst_progress('text', 'Saving source file...');
    bst_progress('inc', 1);
    % Create structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = sqrt(abs(source_diff_dics.pow));
    ResultsMat.nComponents   = 1;
    ResultsMat.Comment       = ['ft_sourceanalysis: ' Method, '_',num2str(f),'Hz'];
    ResultsMat.Function      = Method;
    ResultsMat.Time          = 1;
    ResultsMat.DataFile      = DataFile;
    ResultsMat.HeadModelFile = HeadModelFile;
    ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
    ResultsMat.ChannelFlag   = DataMat.ChannelFlag;
    ResultsMat.GoodChannel   = iChannelsData;
    ResultsMat.SurfaceFile   = HeadModelMat.SurfaceFile;
    ResultsMat.nAvg          = DataMat.nAvg;
    ResultsMat.Leff          = DataMat.Leff;
    ResultsMat.cfg           = source_diff_dics.cfg;
    switch lower(ResultsMat.HeadModelType)
        case 'volume'
            ResultsMat.GridLoc    = HeadModelMat.GridLoc;
            % ResultsMat.GridOrient = [];
        case 'surface'
            ResultsMat.GridLoc    = [];
            % ResultsMat.GridOrient = [];
        case 'mixed'
            ResultsMat.GridLoc    = HeadModelMat.GridLoc;
            ResultsMat.GridOrient = HeadModelMat.GridOrient;
    end
    ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_sourceanalysis: ' Method ' ' Modality]);
    
    % === SAVE OUTPUT FILE ===
    % Output filename
    OutputDir = bst_fileparts(file_fullpath(DataFile));
    ResultFile = bst_process('GetNewFilename', OutputDir, ['results_', Method, '_', Modality, ]);
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    
    % ===== REGISTER NEW FILE =====
    bst_progress('inc', 1);
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = ResultsMat.Comment;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = DataFile;
    newResult.isLink        = 0;
    newResult.HeadModelType = ResultsMat.HeadModelType;
    % Get output study
    sStudyData = bst_get('Study', iStudyData);
    % Add new entry to the database
    iResult = length(sStudyData.Result) + 1;
    sStudyData.Result(iResult) = newResult;
    % Update Brainstorm database
    bst_set('Study', iStudyData, sStudyData);
    % Store output filename
    OutputFiles{end+1} = newResult.FileName;
    % Expand data node
    panel_protocols('SelectNode', [], newResult.FileName);
    %         end
end
% Save database
db_save();
% Hide progress bar
bst_progress('stop');
end



