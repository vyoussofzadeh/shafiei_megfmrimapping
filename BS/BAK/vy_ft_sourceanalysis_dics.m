function varargout = vy_ft_sourceanalysis_dics( varargin )
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
    Index = strfind(HeadModelMat.SurfaceFile, '/');
    subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
    indir = '/MEG_data/epilepsy';
    bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
    bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
    bsdir = fullfile(indir,subj,'brainstorm_db');
    cd(bsdir)
        
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
%     [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
    
    %%
    cfg_main = [];
    cfg_main.outputdir = [];
    cfg_main.freq_of_interest = freq_of_interest;
    cfg_main.sens = datain.grad;
    Run_fft_4dics
    
    %%
    sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
        
    %% Leadfield
    if exist(fullfile(savepath,'leadfield.mat'),'file')~= 2      
        [individual_headmodel, individual_grid] = vy_bs2ft_headmodel(ep_data.all, sourcemodel);
        save(fullfile(savepath,'leadfield.mat'),'individual_headmodel', 'individual_grid');
    else
        load(fullfile(savepath,'leadfield.mat'));
    end
           
    %%
    data = f_data;
    cfg = [];
    cfg.method = 'dics';
    cfg.dics.lambda = '1%';
    cfg.sourcemodel  = individual_grid;
    cfg.frequency    = data.app.freq;
    cfg.headmodel = individual_headmodel;
    cfg.dics.keepfilter = 'yes';
    cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    sourceavg = ft_sourceanalysis(cfg, data.app);
    
    cfg = [];
    cfg.method = 'dics';
    cfg.dics.lambda = '1%';
    cfg.sourcemodel        = individual_grid;
    cfg.sourcemodel.filter = sourceavg.avg.filter;
    cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
    cfg.headmodel = individual_headmodel;
    s_data.bsl      = ft_sourceanalysis(cfg, data.bsl);
    s_data.pst      = ft_sourceanalysis(cfg, data.pst);
    
    %
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
    source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
    source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
    source_diff_dics.pow(source_diff_dics.pow>0)=0;
    source_diff_dics.pow = abs(source_diff_dics.pow);
    
    
    %%
    figure
    m = source_diff_dics.pow;
    bnd.pnt = individual_grid.pos;
    bnd.tri = individual_grid.tri;
    ft_plot_mesh(bnd, 'vertexcolor', sqrt(m));
    colorbar   
   
    %%
    
    % === CREATE OUTPUT STRUCTURE ===
    bst_progress('text', 'Saving source file...');
    bst_progress('inc', 1);
    % Create structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = sqrt(source_diff_dics.pow);
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



