function varargout = process_ft_sourceanalysis_dics( varargin )
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
%% Adding ft
% Command summary goes here
ft_path = '/usr/local/MATLAB_Tools/fieldtrip_20190419';
addpath(ft_path);
ft_defaults


eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% ===== PROCESS =====
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis modified Vahab';
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

h = get(0,'Children');
close(h)

%% Adding extenal functions
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions'))
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/run'))

%%

% ===== GET OPTIONS =====
% Inverse options
Method   = sProcess.options.method.Value{1};
Modality = sProcess.options.sensortype.Value{1};
% Get unique channel files
AllChannelFiles = unique({sInputs.ChannelFile});
% Progress bar
bst_progress('start', 'ft_sourceanalysis', 'Loading input files...', 0, 2*length(sInputs));

% ===== LOOP ON FOLDERS =====
for iChanFile = 1:1%length(AllChannelFiles)
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
    
    % ===== LOOP ON DATA FILES =====
    % Get data files for this channel file
    iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
    % Loop on data files
    for iInput = 1:1%length(iChanInputs)
        
        % === LOAD DATA ===
        % Load data
        DataFile = sInputs(iChanInputs(iInput)).FileName;
        DataMat = in_bst_data(DataFile);
        iStudyData = sInputs(iChanInputs(iInput)).iStudy;
        % Remove bad channels
        iBadChan = find(DataMat.ChannelFlag == -1);
        iChannelsData = setdiff(iChannels, iBadChan);
        % Error: All channels tagged as bad
        if isempty(iChannelsData)
            bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
            return;
        end
        % Convert data file to FieldTrip format
        ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
        % Add data covariance
        ftData.cov = NoiseCovMat.NoiseCov(iChannelsData,iChannelsData);
        % Convert head model to FieldTrip format
        [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
        
        % === CALLING FIELDTRIP FUNCTION ===
        bst_progress('text', 'Calling FieldTrip function: ft_sourceanalysis...');
        % Prepare FieldTrip cfg structure
        cfg           = [];
        cfg.method    = Method;
        cfg.grid      = ftLeadfield;
        cfg.headmodel = ftHeadmodel;
        % Additional options for the method
        switch (Method)
            case 'mne'
                cfg.mne.prewhiten = 'yes';
                cfg.mne.lambda    = 3;
                cfg.mne.scalesourcecov = 'yes';
                Time = DataMat.Time;
                
            case 'lcmv'
                Time = [DataMat.Time(1), DataMat.Time(2)];
                
            case 'dics'
                %% Step1, initial settings
                iChanFile = 1;
                ChannelMat = in_bst_channel(AllChannelFiles{iChanFile});
                iChannels = channel_find(ChannelMat.Channel, Modality);
                iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
                
                %% Step2, reading trials
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
                    trl{iInput} = DataMat.F(iChannels,:);
                    timee {iInput} = DataMat.Time;
                end
                
                ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
                ftData1 = [];
                ftData1.label = ftData.label;
                ftData1.grad = ftData.grad;
                ftData1.dimord = ftData.dimord;
                ftData1.trial = trl;
                ftData1.time = timee;
                
                %% step3, reading headmodel, leadfield, ..
                [sStudyChan, iStudyChan] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
                % Load head model
                HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
                HeadModelMat = in_bst_headmodel(HeadModelFile);
                % Load data covariance matrix
                NoiseCovFile = sStudyChan.NoiseCov(1).FileName;
                NoiseCovMat = load(file_fullpath(NoiseCovFile));
                Index = strfind(HeadModelMat.SurfaceFile, '/');
                subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
                indir = '/MEG_data/epilepsy';
                bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
                bsanatdir = fullfile(indir,subj,'brainstorm_db/anat');
                bsdir = fullfile(indir,subj,'brainstorm_db');
                cd(bsdir)
                
                %% step4, saving path
                Index = strfind(DataFile, '/');
                saveid = DataFile(Index(end)+1:end-4);
                savepath = fullfile(['dics_',saveid]);
                if exist(savepath, 'file') == 0, mkdir(savepath), end
                
                %% step 5: Selecting freq of interest
                disp('Enter the highest freq in the data, e.g. 40')
                fmax = input(':');
                if exist(fullfile(savepath,'tfr.mat'),'file')~= 2
                    % do tfr-decomposition
                    cfg = [];
                    cfg.output     = 'pow';
                    cfg.channel    = 'all';
                    cfg.method     = 'mtmconvol';
                    cfg.taper      = 'hanning';
                    % cfg.taper      = 'dpss';
                    cfg.foi        = 1:2:fmax;
                    cfg.keeptrials = 'yes';
                    cfg.t_ftimwin  = 3./cfg.foi;
                    % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
                    cfg.tapsmofrq  = 0.8 *cfg.foi;
                    cfg.toi        = ftData1.time{1}(1):0.05:ftData1.time{1}(end);
                    tfr        = ft_freqanalysis(cfg, ftData1);
                    save(fullfile(savepath,'tfr.mat'),'tfr');
                else
                    load(fullfile(savepath,'tfr.mat'));
                end
                
                cfg = [];
                cfg.savepath = [];
                cfg.savefile = [];
                cfg.fmax = fmax;
                cfg.toi = [ftData1.time{1}(1), ftData1.time{1}(end)];
                [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
                disp(['time_of_interest:',num2str(time_of_interest),'sec']);
                disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
                L = 0.3;
                
                %% Step 6: selecting time of interest
                datain = ftData1;
                
                toi(1,:) = [-0.3,0];
                toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
                
                disp('===========================================');
                disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),';',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec was selected as contrast intervals']);
                %         disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as post-stim']);
                %         disp('===========================================');
                disp('OK to proceed: 1, No, another time interval: 2');
                ask_time = input(':');
                if ask_time == 2
                    disp('eneter new time interval in sec, e.g., [-0.3,0; 0.7,1.2]');
                    toi = input(':');
                end
                
                %
                ep_data = vy_epoch(datain, toi);
                cfg = [];
                ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
                
                %% step7: grand average
                %             cfg = [];
                %             cfg.layout = 'neuromag306mag.lay';
                %             lay = ft_prepare_layout(cfg);
                %             Run_grandmean
                
                %% step8: spectral analysis
                cfg_main = [];
                cfg_main.sens = datain.grad;
                cfg_main.outputdir = savepath;
                cfg_main.freq_of_interest  = freq_of_interest; % Hz
                Run_fft_4dics
                
                %% step8: Surface-based source analysis, from BS
                sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
                [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
                
                %         figure;
                %         % ft_plot_vol(individual_headmodel_surf, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
                %         %         ft_plot_vol(ftHeadmodel , 'facecolor', 'cortex', 'edgecolor', 'none');     alpha 0.5;
                %         hold on;
                %         ft_plot_headshape(sourcemodel);
                %         ft_plot_mesh(ftLeadfield.pos(ftLeadfield.inside, :));
                %         view ([-10 40 10])
                %
                cfg = [];
                cfg.method = 'dics';
                cfg.dics.lambda = '1%';
                cfg.sourcemodel  = ftLeadfield;
                cfg.frequency    = f_data.app.freq;
                cfg.headmodel = ftHeadmodel;
                cfg.dics.keepfilter = 'yes';
                cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
                sourceavg = ft_sourceanalysis(cfg, f_data.app);
                
                cfg = [];
                cfg.method = 'dics';
                cfg.dics.lambda = '1%';
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
                source_diff_dics.pow = abs(source_diff_dics.pow);
                
                %         figure
                %         m = source_diff_dics.pow;
                %         bnd.pnt = sourcemodel.pos;
                %         bnd.tri = sourcemodel.tri;
                %         ft_plot_mesh(bnd, 'vertexcolor', abs(m));
                %         colorbar
                %
                
            case 'pcc'
                % % ft_freqanalysis %
                % cfg=[];
                % cfg.method      = 'mtmfft';
                % cfg.output      = 'fourier';  % gives the complex Fourier spectra
                % cfg.foilim      = [60 60];    % analyse 40-80 Hz (60 Hz +/- 20 Hz smoothing)
                % cfg.taper       = 'dpss';
                % cfg.tapsmofrq   = 20;
                % cfg.keeptrials  = 'yes';      % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
                % cfg.keeptapers  = 'yes';
                % freq = ft_freqanalysis(cfg, data);
                %
                % % compute common spatial filter AND project all trials through it %
                % cfg=[];
                % cfg.method      = 'pcc';
                % cfg.grid        = grid;       % previously computed grid
                % cfg.headmodel   = vol;        % previously computed volume conduction model
                % cfg.frequency   = 60;
                % cfg.keeptrials  = 'yes';      % keep single trials. Only necessary if you are interested in reconstructing single trial data
                % source = ft_sourceanalysis(cfg, freq);
        end
        % Call FieldTrip function
        %             ftSource = ft_sourceanalysis(cfg, ftData);
        
        % === CREATE OUTPUT STRUCTURE ===
        bst_progress('text', 'Saving source file...');
        bst_progress('inc', 1);
        % Create structure
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        ResultsMat.ImageGridAmp  = sqrt((source_diff_dics.pow));
        ResultsMat.nComponents   = 1;
        ResultsMat.Comment       = ['ft_sourceanalysis: ' Method, '_',num2str(f),'Hz_', num2str(toi(2,1)),'_',num2str(toi(2,2)), 'sec ', datestr(now, 'dd/mm/yy-HH:MM')];
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
    end
end
% Save database
db_save();
% Hide progress bar
bst_progress('stop');
end