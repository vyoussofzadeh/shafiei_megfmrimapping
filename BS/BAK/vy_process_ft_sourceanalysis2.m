function varargout = vy_process_ft_sourceanalysis( varargin )
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
sProcess.Index       = 356;
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
    
    %% tfr analysis
    
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
    
    cfg = [];
    cfg.savepath = [];
    cfg.savefile = [];
    [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
    disp(['time_of_interest:',num2str(time_of_interest),'sec']);
    disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
    L = 0.3;
    
    %% Selecting time of interest
    datain = ftData1;
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
    % Run_time
    ep_data = vy_epoch(datain, toi);
    
    %- Appending data
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %%
    [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
    
    %%
    cfg_main = [];
    cfg_main.outputdir = [];
    cfg_main.freq_of_interest = freq_of_interest;
    cfg_main.sens = ftData1.grad;
    Run_fft_4dics
    
    %%
    Index = strfind(HeadModelMat.SurfaceFile, '/');
    subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
    
    indir = '/MEG_data/epilepsy';
    bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
    bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
    
    sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
    
    %%
%     anatfile = fullfile(bsanatdir, subj);
%     
%     disp(anatfile)
%     d = rdir(fullfile(anatfile,'subjectimage*.mat'));
%     % d = rdir(fullfile(datafile,'*ras.mat'));
%     if ~isempty(d)
%         sMRI = d.name;
%         %     cd(anatfile)
%         %     cd ..
%         OutputFile = 'T1.nii';
%         out_mri_nii(sMRI, OutputFile);
%     end
    
    %%
%     [individual_headmodel, individual_grid] = vy_bs2ft_headmodel(ep_data.all, sourcemodel);
    % save('test.mat','individual_headmodel', 'individual_grid');
    load('test.mat');
    
    %% Volumetric-based analysis   
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
%     
%     %
%     cfg = [];
%     cfg.method = 'singleshell';
%     headmodel = ft_prepare_headmodel(cfg, outanat.individual_seg);
%     figure,ft_plot_headmodel(headmodel);
%     headmodel = ft_convert_units(headmodel, 'm');
%     
% %     headmodel1 = ft_transform_geometry(ChannelMat.TransfMeg{1}./outanat.individual_seg.transform,headmodel);
%     headmodel1 = ft_transform_geometry(ChannelMat.TransfMeg{2},headmodel);
%     figure,ft_plot_headmodel(headmodel1);
%     
%     %
%     individual_grid = ft_convert_units(individual_grid, 'm');
%     individual_headmodel = ft_convert_units(individual_headmodel, 'm');
% 
%     figure;
%     ft_plot_headmodel(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     hold on;
%     ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
%     view ([-10 40 10])
%     
%     a = [];
%     a.pos = ChannelMat.HeadPoints.Loc'; % from channle vectorview
%     a.label = ChannelMat.HeadPoints.Label;
%     
%     figure; hold on;
% %     ft_plot_headmodel(ftHeadmodel);
%     ft_plot_headmodel(headmodel1);
% %     ft_plot_headmodel(individual_headmodel);
%     alpha 0.5;
%     ft_plot_mesh(individual_grid, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
%     hold on;
%     ft_plot_headshape(sourcemodel.pos);

    
    %%
    data = f_data;
    cfg = [];
    cfg.method = 'dics';
    cfg.dics.lambda = '1%';
    cfg.sourcemodel           = individual_grid;
    cfg.frequency    = data.app.freq;
    % cfg.sourcemodel = ftLeadfield;
    % cfg.grid = cfg_main.grid;
    cfg.headmodel = individual_headmodel;
%     cfg.headmodel = headmodel1;
%     cfg.headmodel = ftHeadmodel;
    cfg.dics.keepfilter = 'yes';
    cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    sourceavg = ft_sourceanalysis(cfg, data.app);
    
    cfg = [];
    cfg.method = 'dics';
    cfg.dics.lambda = '1%';
    cfg.sourcemodel        = individual_grid;
    % cfg.sourcemodel = ftLeadfield;
    % cfg.grid = cfg_main.grid;
    cfg.sourcemodel.filter = sourceavg.avg.filter;
    cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    cfg.headmodel = individual_headmodel;
%     cfg.headmodel = headmodel1;
%     cfg.headmodel = ftHeadmodel;
    s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
    % s_data.bsl.avg.pow = s_data.bsl.avg.pow./max(s_data.bsl.avg.pow(:));
    s_data.pst      = ft_sourceanalysis(cfg, data.pst);
    % s_data.pst.avg.pow = s_data.pst.avg.pow./max(s_data.pst.avg.pow(:));
    
    %
    cfg = [];
    cfg.parameter = 'pow';
%     cfg.operation = '(x1-x2)/(x1+x2)';
    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
    source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
    % source_diff_dics.pos     = cfg_main.template_grid.pos;
    % source_diff_dics.dim     = cfg_main.template_grid.dim;
    % source_diff_dics.inside  = cfg_main.template_grid.inside;
    % source_diff_dics.pow(source_diff_dics.pow>0)=0;
    source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
    source_diff_dics.pow(source_diff_dics.pow>0)=0;
    source_diff_dics.pow = abs(source_diff_dics.pow);
    
    
%     %%
%     figure
%     m = source_diff_dics.pow;
%     bnd.pnt = individual_grid.pos;
%     bnd.tri = individual_grid.tri;
%     ft_plot_mesh(bnd, 'vertexcolor', sqrt(m));
%     colorbar   
    
    %%
    % data = f_data;
    % cfg = [];
    % cfg.method = 'dics';
    % cfg.dics.lambda = '5%';
    % cfg.frequency    = data.app.freq;
    % cfg.sourcemodel = ftLeadfield;
    % % cfg.grid = cfg_main.grid;
    % cfg.headmodel = ftHeadmodel;
    % cfg.dics.keepfilter = 'yes';
    % cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    % sourceavg = ft_sourceanalysis(cfg, data.app);
    %
    % cfg = [];
    % cfg.method = 'dics';
    % cfg.dics.lambda = '5%';
    % cfg.sourcemodel = ftLeadfield;
    % % cfg.grid = cfg_main.grid;
    % cfg.sourcemodel.filter = sourceavg.avg.filter;
    % cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    % cfg.headmodel = ftHeadmodel;
    % s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
    % % s_data.bsl.avg.pow = s_data.bsl.avg.pow./max(s_data.bsl.avg.pow(:));
    % s_data.pst      = ft_sourceanalysis(cfg, data.pst);
    % % s_data.pst.avg.pow = s_data.pst.avg.pow./max(s_data.pst.avg.pow(:));
    %
    % %%
    % cfg = [];
    % cfg.parameter = 'pow';
    % % cfg.operation = '(x1-x2)/(x1+x2)';
    % cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
    % source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
    % % source_diff_dics.pos     = cfg_main.template_grid.pos;
    % % source_diff_dics.dim     = cfg_main.template_grid.dim;
    % % source_diff_dics.inside  = cfg_main.template_grid.inside;
    % % source_diff_dics.pow(source_diff_dics.pow>0)=0;
    % source_diff_dics.pow(source_diff_dics.pow>0)=0;
    % source_diff_dics.pow = abs(source_diff_dics.pow);
    %
    % %%
    % figure
    % m = source_diff_dics.pow;
    % bnd.pnt = ftLeadfield.pos;
    % bnd.tri = sourcemodel.tri;
    % ft_plot_mesh(bnd, 'vertexcolor', m);
    % colorbar
    
    
    % cfg = [];
    % cfg.mask = 'pow';
    % cfg.loc = 'min';
    % cfg.template = cfg_main.template_mri;
    % cfg.savefile = savefig;
    % cfg.volnorm     = 2; % yes: 1
    % source_int_dics = vy_source_plot(cfg, source_diff_dics);
    
    %%
    % source_diff_dics.tri = sourcemodel.tri;
    % cfg = [];
    % cfg.method          = 'surface';
    % cfg.funparameter    = 'pow';
    % cfg.colorbar        = 'no';
    % cfg.opacitymap    = 'rampdown';
    % cfg.funcolormap        = brewermap(256, '*RdYlBu');
    % ft_sourceplot(cfg, source_diff_dics);
    % colorbar
    
    
    %%
    % cfg = [];
    % cfg.headmodel = ftHeadmodel;
    % cfg.sourcemodel = ftLeadfield;
    % % cfg.grid = cfg_main.grid;
    % cfg.mtag = 'dics';
    % s_data_dics = vy_source_freq(cfg, f_data);
    
    
    %%
    
    % ===== LOOP ON DATA FILES =====
    % Get data files for this channel file
    %         iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
    %         % Loop on data files
    %         for iInput = 1:length(iChanInputs)
    %
    %             % === LOAD DATA ===
    %             % Load data
    %             DataFile = sInputs(iChanInputs(iInput)).FileName;
    %             DataMat = in_bst_data(DataFile);
    %             iStudyData = sInputs(iChanInputs(iInput)).iStudy;
    %             % Remove bad channels
    %             iBadChan = find(DataMat.ChannelFlag == -1);
    %             iChannelsData = setdiff(iChannels, iBadChan);
    %             % Error: All channels tagged as bad
    %             if isempty(iChannelsData)
    %                 bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
    %                 return;
    %             end
    %             % Convert data file to FieldTrip format
    %             ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
    %             % Add data covariance
    %             ftData.cov = NoiseCovMat.NoiseCov(iChannelsData,iChannelsData);
    %             % Convert head model to FieldTrip format
    %             [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
    %
    %             % === CALLING FIELDTRIP FUNCTION ===
    %             bst_progress('text', 'Calling FieldTrip function: ft_sourceanalysis...');
    %             % Prepare FieldTrip cfg structure
    %             cfg           = [];
    %             cfg.method    = Method;
    % %             cfg.grid      = ftLeadfield;
    %             cfg.sourcemodel      = ftLeadfield;
    %             cfg.headmodel = ftHeadmodel;
    %             % Additional options for the method
    %
    %
    %             %%
    %             figure; hold on;
    % % ft_plot_vol(ftHeadmodel); alpha 0.5;
    % ft_plot_mesh(ftLeadfield, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
    % % ft_plot_headshape(sourcemodel);
    %
    % a = [];
    % a.pos = ChannelMat.HeadPoints.Loc';
    % a.label = ChannelMat.HeadPoints.Label;
    % ft_plot_headshape(a);
    
    
    % ft_plot_sens(ftData.grad, 'unit', 'm', 'coilsize', 10, 'chantype', 'meggrad');
    % grad = ft_convert_units(ep_data1.pst.grad,'m');ft_plot_sens(ChannelMat.HeadPoints)
    
    
    %             %%
    %             a = Method(1);
    %             switch a{1}
    %                 case 'mne'
    %                     cfg = [];
    %                     cfg.method = 'mne';
    %                     cfg.mne.prewhiten = 'yes';
    %                     cfg.sourcemodel      = ftLeadfield;
    %                     cfg.headmodel = ftHeadmodel;
    %                     cfg.mne.lambda    = 3;
    %                     cfg.mne.scalesourcecov = 'yes';
    %                     Time = DataMat.Time;
    %
    %                 case 'lcmv'
    %
    %                     cfg.method      = 'lcmv';
    %                     Time = [DataMat.Time(1), DataMat.Time(2)];
    %
    %                 case 'dics'
    %
    %                     cfg = [];
    %                     cfg.method    = 'mtmfft';
    % %                     cfg.output    = 'powandcsd';
    %                     cfg.output       = 'fourier';
    %                     cfg.tapsmofrq = 4;
    %                     cfg.foilim    = [10 10];
    %                     freqPre = ft_freqanalysis(cfg, ftData);
    %
    %                     cfg              = [];
    %                     cfg.method       = 'mtmfft';
    %                     % cfg.output       = 'pow';
    %                     % cfg.pad          = 'nextpow2';
    %                     cfg.output       = 'fourier';
    %                     cfg.keeptrials   = 'yes';
    %                     cfg.foilim = [2 40];
    % %                     cfg.foilim    = [18 18];
    %                     cfg.taper    = 'hanning';
    %                     cfg.taper    = 'dpss';
    %                     cfg.tapsmofrq    = 4;
    %                     cfg.pad          = 4;
    %                     freqPre1             = ft_freqanalysis(cfg, ftData);
    %
    %                     psd = squeeze(mean(mean(abs(freqPre1.fourierspctrm),2),1));
    %                     ff = linspace(1, cfg.foilim(2), length(psd));
    %                     figure,plot(ff,psd)
    %                     xlabel('Hz'); ylabel('psd')
    %
    %                     % EXAMPLE 1
    %                     cfg                = [];
    %                     cfg.grid           = ftLeadfield;
    %                     cfg.frequency      = freqPre.freq;
    %                     cfg.vol            = ftHeadmodel;
    %                     cfg.gradfile       = freqPre.grad;
    %                     cfg.projectnoise   = 'yes';
    %                     cfg.keeptrials     = 'no';
    %                     cfg.keepfilter     = 'yes';
    %                     cfg.keepcsd        = 'yes';
    %                     cfg.keepmom        = 'yes';
    %                     cfg.lambda         = '5%';
    % %                     cfg.lambda         = 0.1 * mean(f.powspctrm(:,nearest(cfg.frequency)),1);
    %                     cfg.method         = 'dics';
    %                     cfg.feedback       = 'textbar';
    %                     ftSource             = ft_sourceanalysis(cfg,freqPre);
    %
    %                     % EXAMPLE 2
    %                     % % freqanalysis %
    %                     % cfg=[];
    %                     % cfg.method      = 'mtmfft';
    %                     % cfg.output      = 'powandcsd';  % gives power and cross-spectral density matrices
    %                     % cfg.foilim      = [60 60];      % analyse 40-80 Hz (60 Hz +/- 20 Hz smoothing)
    %                     % cfg.taper       = 'dpss';
    %                     % cfg.tapsmofrq   = 20;
    %                     % cfg.keeptrials  = 'yes';        % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
    %                     % cfg.keeptapers  = 'no';
    %                     %
    %                     % freq = ft_freqanalysis(cfg, data);
    %                     %
    %                     % % compute common spatial filter %
    %                     % cfg=[];
    %                     % cfg.method      = 'dics';
    %                     % cfg.grid        = grid;         % previously computed grid
    %                     % cfg.headmodel   = vol;          % previously computed volume conduction model
    %                     % cfg.frequency   = 60;
    %                     % cfg.dics.keepfilter  = 'yes';        % remember the filter
    %                     %
    %                     % source = ft_sourceanalysis(cfg, freq);
    %
    %                 case 'pcc'
    %                     % % ft_freqanalysis %
    %                     % cfg=[];
    %                     % cfg.method      = 'mtmfft';
    %                     % cfg.output      = 'fourier';  % gives the complex Fourier spectra
    %                     % cfg.foilim      = [60 60];    % analyse 40-80 Hz (60 Hz +/- 20 Hz smoothing)
    %                     % cfg.taper       = 'dpss';
    %                     % cfg.tapsmofrq   = 20;
    %                     % cfg.keeptrials  = 'yes';      % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
    %                     % cfg.keeptapers  = 'yes';
    %                     % freq = ft_freqanalysis(cfg, data);
    %                     %
    %                     % % compute common spatial filter AND project all trials through it %
    %                     % cfg=[];
    %                     % cfg.method      = 'pcc';
    %                     % cfg.grid        = grid;       % previously computed grid
    %                     % cfg.headmodel   = vol;        % previously computed volume conduction model
    %                     % cfg.frequency   = 60;
    %                     % cfg.keeptrials  = 'yes';      % keep single trials. Only necessary if you are interested in reconstructing single trial data
    %                     % source = ft_sourceanalysis(cfg, freq);
    %             end
    %             % Call FieldTrip function
    %             %             ftSource = ft_sourceanalysis(cfg, ftData);
    %
    %             %%
    %
    %             Index = strfind(HeadModelMat.SurfaceFile, '/');
    %             subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
    %
    %             indir = '/MEG_data/epilepsy';
    %             bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
    %             bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
    %
    %             sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
    
    %             figure;
    %             bnd.pnt = sourcemodel.pos;
    %             bnd.tri = sourcemodel.tri;
    %             %             bnd.funcolormap =  brewermap(256, '*RdYlBu');
    %             cfg.opacitymap    = 'rampdown';
    %             ft_plot_mesh(bnd, 'vertexcolor', sqrt(ftSource.avg.pow));
    %             colorbar
    %             light ('Position',[-70 20 50])
    %             material dull
    
    
    %             ftSource.tri = sourcemodel.tri;
    %             cfg = [];
    %             cfg.method          = 'surface';
    %             cfg.funparameter    = 'pow';
    %             cfg.colorbar        = 'no';
    %             cfg.opacitymap    = 'rampdown';
    %             cfg.funcolormap        = brewermap(256, '*RdYlBu');
    %             ft_sourceplot(cfg, ftSource);
    
    %%
    
    % === CREATE OUTPUT STRUCTURE ===
    bst_progress('text', 'Saving source file...');
    bst_progress('inc', 1);
    % Create structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = sqrt(source_diff_dics.pow);
    ResultsMat.nComponents   = 1;
    ResultsMat.Comment       = ['ft_sourceanalysis: ' Method];
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



