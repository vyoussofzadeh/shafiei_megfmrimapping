function fcn_hcp_meg_process_connectivity(hcp_dir, subjList, badChannels, reports_dir)
% hcp_meg_process: Script to process Resting-state MEG data from Human
% Connectime Project
%
% Modified from Brainstorm online tutorials:
%     https://neuroimage.usc.edu/brainstorm/Tutorials/HCP-MEG
%
% Inputs:
%     hcp_dir: Directory with subject subdirectories that have the unzipped
%     HCP files
%     subjList: List of subjects in a Matlab cell format (e.g. {'102816'})
%     badChannels: List of bad channels for each subject in a Matlab cell
%     fromat (e.g. {{'A2', 'A237', 'A244', 'A246', 'A8'}, {'A126', 'A2',
%     'A244', 'A246'}})
%     reports_dir: Directory to save Brainstorm generated reports for each
%     subject
%
% Outputs:
%     saves vertex-level, band-passed filtered AEC and PLV connectivity
%     matrices for each subject
%
% @=============================================================================
% This code is a part of meg-fmri mapping research project and it is
% still under development.
% =============================================================================@
%
% Author: Golia Shafiei, 2021

%% ===== SCRIPT VARIABLES =====
% Full list of subjects to process, bad channels, frequency bands
SubjectNames = subjList;
BadChannels = badChannels;

freqL.delta = 2.;
freqH.delta = 4.;
freqL.theta = 5.;
freqH.theta = 7.;
freqL.alpha = 8.;
freqH.alpha = 12.;
freqL.beta = 15.;
freqH.beta = 29.;
freqL.lgamma = 30.;
freqH.lgamma = 59.;
freqL.hgamma = 60.;
freqH.hgamma = 90.;

%% ===== CREATE PROTOCOL =====
% The protocol name has to be a valid folder name (no spaces, no weird characters...)
ProtocolName = 'megHCP';
% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
    %brainstorm server
end
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

%% ===== FILES TO IMPORT =====
% You have to specify the folder in which the tutorial dataset is unzipped
if (nargin < 3) || isempty(hcp_dir) || ~file_exist(hcp_dir)
    error('The first argument must be the full path to the tutorial dataset folder.');
end
% Output folder for reports
if (nargin < 3) || isempty(reports_dir) || ~isdir(reports_dir)
    reports_dir = [];
end

%% ===== PRE-PROCESS AND IMPORT =====
for iSubj = 1:length(SubjectNames)
    tic
    % Start a new report (one report per subject)
    bst_report('Start');
    fprintf('\n===== IMPORT: SUBJECT #%d =====\n', iSubj);

    % If subject already exists: delete it
    [sSubject, iSubject] = bst_get('Subject', SubjectNames{iSubj});
    if ~isempty(sSubject)
        db_delete_subjects(iSubject);
    end

    % ===== FILES TO IMPORT =====
    % Build the path of the files to import
    AnatDir    = fullfile(hcp_dir, SubjectNames{iSubj}, 'MEG', 'anatomy');
    Run1File   = fullfile(hcp_dir, SubjectNames{iSubj}, 'unprocessed', 'MEG', '3-Restin', '4D', 'c,rfDC');
    NoiseFile  = fullfile(hcp_dir, SubjectNames{iSubj}, 'unprocessed', 'MEG', '1-Rnoise', '4D', 'c,rfDC');
    % Check if the folder contains the required files
    if ~file_exist(AnatDir) || ~file_exist(Run1File) || ~file_exist(NoiseFile)
        error(['The folder ' hcp_dir ' does not contain subject #' SubjectNames{iSubj} 'from the HCP-MEG distribution.']);
    end

    % ===== IMPORT DATA =====
    % Process: Import anatomy folder
    bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectNames{iSubj}, ...
        'mrifile',     {AnatDir, 'HCPv3'}, ...
        'nvertices',   15000);

    % Process: Create link to raw files
    sFilesRun1 = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',  SubjectNames{iSubj}, ...
        'datafile',     {Run1File, '4D'}, ...
        'channelalign', 1);
    sFilesNoise = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',  SubjectNames{iSubj}, ...
        'datafile',     {NoiseFile, '4D'}, ...
        'channelalign', 1);
    sFilesRaw = [sFilesRun1, sFilesNoise];

    % Process: Resample: 508.63Hz
    sFilesResamp = bst_process('CallProcess', 'process_resample', sFilesRaw, [], ...
        'freq',     508.6275, ...
        'read_all', 1);

    % ===== PRE-PROCESSING =====
    % Process: Notch filter: 60Hz 120Hz 180Hz 240Hz 300Hz
    sFilesNotch = bst_process('CallProcess', 'process_notch', sFilesResamp, [], ... %sFilesRaw
        'freqlist',    [60, 120, 180, 240, 300], ...
        'sensortypes', 'MEG, EEG', ...
        'read_all',    1);

    % Process: High-pass:0.3Hz
    sFilesBand = bst_process('CallProcess', 'process_bandpass', sFilesNotch, [], ...
        'sensortypes', 'MEG, EEG', ...
        'highpass',    0.3, ...
        'lowpass',     0, ...
        'attenuation', 'strict', ...  % 60dB
        'mirror',      0, ...
        'useold',      0, ...
        'read_all',    1);

    % Process: Power spectrum density (Welch)
    sFilesPsdAfter = bst_process('CallProcess', 'process_psd', sFilesBand, [], ...
        'timewindow',  [], ...
        'win_length',  4, ...
        'win_overlap', 50, ...
        'sensortypes', 'MEG, EEG', ...
        'edit',        struct(...
             'Comment',         'Power', ...
             'TimeBands',       [], ...
             'Freqs',           [], ...
             'ClusterFuncTime', 'none', ...
             'Measure',         'power', ...
             'Output',          'all', ...
             'SaveKernel',      0));

    % Mark bad channels
    bst_process('CallProcess', 'process_channel_setbad', sFilesBand, [], ...
                'sensortypes', BadChannels{iSubj});

    % Process: Snapshot: Frequency spectrum
    bst_process('CallProcess', 'process_snapshot', sFilesPsdAfter, [], ...
        'target',         10, ...  % Frequency spectrum
        'modality',       1);      % MEG (All)

    % Process: Delete folders
    bst_process('CallProcess', 'process_delete', ...
        [sFilesRaw, sFilesNotch, sFilesResamp], [], ...
        'target', 2);  % Delete folders

    % ===== ARTIFACT CLEANING =====
    % Process: Select data files in: */*
    sFilesBand = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname', SubjectNames{iSubj});

    % Process: Select file names with tag: 3-Restin
    sFilesRest = bst_process('CallProcess', 'process_select_tag', sFilesBand, [], ...
        'tag',    '3-Restin', ...
        'search', 1, ...  % Search the file names
        'select', 1);  % Select only the files with the tag

    % Process: Detect heartbeats
    bst_process('CallProcess', 'process_evt_detect_ecg', sFilesRest, [], ...
        'channelname', 'ECG+, -ECG-', ...
        'timewindow',  [], ...
        'eventname',   'cardiac');

    % Process: Detect eye blinks
    bst_process('CallProcess', 'process_evt_detect_eog', sFilesRest, [], ...
        'channelname', 'VEOG+, -VEOG-', ...
        'timewindow', [], ...
        'eventname', 'blink');

    % Process: Remove simultaneous (keep blinks over heart beats)
    bst_process('CallProcess', 'process_evt_remove_simult', sFilesRest, [], ...
        'remove', 'cardiac', ...
        'target', 'blink', ...
        'dt', 0.25, ...
        'rename', 0);

    % Process: SSP ECG: cardiac (force remove 1st component)
    bst_process('CallProcess', 'process_ssp_ecg', sFilesRest, [], ...
        'eventname',   'cardiac', ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'select',      1);

    % Process: SSP EOG: blink (force remove 1st component)
    bst_process('CallProcess', 'process_ssp_eog', sFilesRest, [], ...
        'eventname', 'blink', ...
        'sensortypes', 'MEG', ...
        'usessp', 1, ...
        'select', 1);

    % SSP: Noisy signal, Sacades, EMG
    % Process: Detect other artifacts (mark noisy segments)
    bst_process('CallProcess', 'process_evt_detect_badsegment', ...
        sFilesRest, [], ...
        'timewindow', [], ...
        'sensortypes', 'MEG, EEG', ...
        'threshold', 3, ...  % 3
        'isLowFreq', 1, ...
        'isHighFreq', 1);

    % Process: SSP for low frequencies (saccades) 1 - 7 Hz (remove 1st)
    bst_process('CallProcess', 'process_ssp', sFilesRest, [], ...
        'timewindow',  [], ...
        'eventname',   '1-7Hz', ...
        'eventtime',   [], ...
        'bandpass',    [1.5, 7], ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);

    % Process: SSP for high frequencies (muscle) 40 - 240 Hz (remove 1st)
    bst_process('CallProcess', 'process_ssp', sFilesRest, [], ...
        'timewindow',  [], ...
        'eventname',   '40-240Hz', ...
        'eventtime',   [], ...
        'bandpass',    [40, 240], ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);

    % Process: Snapshot: Sensors/MRI registration
    bst_process('CallProcess', 'process_snapshot', sFilesRest, [], ...
        'target',         1, ...  % Sensors/MRI registration
        'modality',       1, ...  % MEG (All)
        'orient',         1);  % left

    % Process: Snapshot: SSP projectors
    bst_process('CallProcess', 'process_snapshot', sFilesRest, [], ...
        'target',         2, ...  % SSP projectors
        'modality',       1);     % MEG (All)

    % ===== SOURCE ESTIMATION =====
    % Process: Select file names with tag: task-rest
    sFilesNoise = bst_process('CallProcess', 'process_select_tag', sFilesBand, [], ...
        'tag',    '1-Rnoise', ...
        'search', 1, ...  % Search the file names
        'select', 1);  % Select only the files with the tag

%    To save the full time series uncomment these, and replace sFilesNoise
%    and sFilesRest with sFilesNoiseFull and sFilesRestFull:
    sFilesNoiseFull = bst_process('CallProcess', 'process_import_data_time', sFilesNoise, [], ...
        'subjectname',  SubjectNames{iSubj}, ...
        'condition',    '', ...
        'datafile',     {'', ''}, ...
        'timewindow',   [], ...
        'split',        0, ...
        'ignoreshort',  0, ...
        'channelalign', 0, ...
        'usectfcomp',   0, ...
        'usessp',       0, ...
        'freq',         [], ...
        'baseline',     []);

    sFilesRestFull = bst_process('CallProcess', 'process_import_data_time', sFilesRest, [], ...
        'subjectname',  SubjectNames{iSubj}, ...
        'condition',    '', ...
        'datafile',     {'', ''}, ...
        'timewindow',   [], ...
        'split',        0, ...
        'ignoreshort',  0, ...
        'channelalign', 0, ...
        'usectfcomp',   0, ...
        'usessp',       0, ...
        'freq',         [], ...
        'baseline',     []);

    % Process: Compute covariance (noise or data)
    bst_process('CallProcess', 'process_noisecov',  sFilesNoiseFull, [], ...
        'baseline',       [], ...
        'sensortypes',    'MEG', ...
        'target',         1, ...  % Noise covariance     (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       1, ...
        'copysubj',       0, ...
        'replacefile',    1);  % Replace

    temp = sFilesRest.FileName;
    sTime = load(file_fullpath(temp), 'Time');
    bst_process('CallProcess', 'process_noisecov', sFilesRestFull, [], ...
        'baseline',       [sTime.Time(1) sTime.Time(end)], ...
        'datatimewindow', [sTime.Time(1) sTime.Time(end)], ...
        'sensortypes',    'MEG', ...
        'target',         2, ...  % Data covariance      (covariance over data time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       1, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace

    % Process: Compute head model
    bst_process('CallProcess', 'process_headmodel', sFilesRestFull, [], ...
        'sourcespace', 1, ...  % Cortex surface
        'meg',         3);     % Overlapping spheres

    % Process: Compute sources [2018]: LCMV
    sSrcRest = bst_process('CallProcess', 'process_inverse_2018', sFilesRestFull, [], ...
        'output',  2, ...  % Kernel only: one per file: 2; Full results: 3
        'inverse', struct(...
             'Comment',        'PNAI: MEG', ...
             'InverseMethod',  'lcmv', ...
             'InverseMeasure', 'nai', ...
             'SourceOrient',   {{'fixed'}}, ...
             'Loose',          0.2, ...
             'UseDepth',       1, ...
             'WeightExp',      0.5, ...
             'WeightLimit',    10, ...
             'NoiseMethod',    'median', ...
             'NoiseReg',       0.1, ...
             'SnrMethod',      'rms', ...
             'SnrRms',         1e-06, ...
             'SnrFixed',       3, ...
             'ComputeKernel',  1, ... % change to 1 for Kernel only and to 0 for Full results
             'DataTypes',      {{'MEG'}}));
         
%     % Process: Compute sources [2018]: sLoreta
%     sSrcRest = bst_process('CallProcess', 'process_inverse_2018', sFilesRestFull, [], ...
%         'output',  2, ...  % Kernel only: one per file
%         'inverse', struct(...
%              'Comment',        'sLORETA: MEG', ...
%              'InverseMethod',  'minnorm', ...
%              'InverseMeasure', 'sloreta', ...
%              'SourceOrient',   {{'fixed'}}, ...
%              'Loose',          0.2, ...
%              'UseDepth',       0, ...
%              'WeightExp',      0.5, ...
%              'WeightLimit',    10, ...
%              'NoiseMethod',    'median', ...
%              'NoiseReg',       0.1, ...
%              'SnrMethod',      'fixed', ...
%              'SnrRms',         1e-06, ...
%              'SnrFixed',       3, ...
%              'ComputeKernel',  1, ...
%              'DataTypes',      {{'MEG'}}));
     

     % get resolution matrix for localization error estimation
     sourcedata = in_bst_results(sSrcRest.FileName);
     
     sStudy = bst_get('Study');
     headmodel = sStudy.HeadModel.FileName;
     sHeadModel = in_bst_headmodel(headmodel);
     Gain_constrained = bst_gain_orient(sHeadModel.Gain, sHeadModel.GridOrient);
     gain = Gain_constrained(sourcedata.GoodChannel, :);
     
     resolutionMat = sourcedata.ImagingKernel*gain;
     
     % Save Vertex x Vertex Resolution Matrix
     outpath = fullfile(hcp_dir, 'brainstormResults', ...
         'resolutionMatrix_lcmv', SubjectNames{iSubj});
     if ~isfolder(outpath)
         mkdir(outpath)
         save(fullfile(outpath, strcat(SubjectNames{iSubj}, ...
             '_resMat.mat')), 'resolutionMat')
     else
         save(fullfile(outpath, strcat(SubjectNames{iSubj}, ...
             '_resMat.mat')), 'resolutionMat')
     end

     % ===== Parcellation =====
     % Load Atlas
     ScoutFile = load(fullfile(hcp_dir, 'parcellationData', 'parcellation', ...
         'scout_schaefer400.4k.label_401.mat'));

     % Load Subj Atlas
     sSubject = bst_get('Subject', SubjectNames{iSubj});
     CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
     sCortex = in_tess_bst(CortexFile);
     idxAtlas = length(sCortex.Atlas)+1;

     % Update Subj Atlas
     sCortex.Atlas(idxAtlas).Name = ScoutFile.Name;
     sCortex.Atlas(idxAtlas).Scouts = ScoutFile.Scouts;

     % Save Updated Atlas
     bst_save(file_fullpath(CortexFile), sCortex, 'v7');
    
     % Save Vertex x 3 coordinates for resolution matrix
     vertices = sCortex.Vertices;
     outpath = fullfile(hcp_dir, 'brainstormResults', ...
         'resolutionMatrix_coordinates', SubjectNames{iSubj});
     if ~isfolder(outpath)
         mkdir(outpath)
         save(fullfile(outpath, strcat(SubjectNames{iSubj}, ...
             '_vertexCoor.mat')), 'vertices')
     else
         save(fullfile(outpath, strcat(SubjectNames{iSubj}, ...
             '_vertexCoor.mat')), 'vertices')
     end

     % ===== Connectivity & BAND-PASS FILTERING: FREQ BANDS =====
     % Frequency bands
     sConnBandAEC = struct('delta', ' ', 'theta', ' ', 'alpha', ' ', ...
         'beta', ' ', 'lgamma', ' ', 'hgamma', ' ');
     sConnBandPLV = struct('delta', ' ', 'theta', ' ', 'alpha', ' ', ...
         'beta', ' ', 'lgamma', ' ', 'hgamma', ' ');
     sBandLabels = fieldnames(sConnBandPLV);

     for iBand = 1:numel(sBandLabels)
         % Process: Connectivity AEC - PCA
         sConnBandAEC.(sBandLabels{iBand}) = bst_process('CallProcess', 'process_aec1n', ...
             sSrcRest, [], ...
             'timewindow', [sTime.Time(1) sTime.Time(end)], ...
             'scouts',     {'schaefer400.4k.label', {'0', '1', '10', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '11', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '12', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '13', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '14', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '15', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '16', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '17', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '18', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '19', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '2', '20', '200', '201', '202', '203', '204', '205', '206', '207', '208', '209', '21', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219', '22', '220', '221', '222', '223', '224', '225', '226', '227', '228', '229', '23', '230', '231', '232', '233', '234', '235', '236', '237', '238', '239', '24', '240', '241', '242', '243', '244', '245', '246', '247', '248', '249', '25', '250', '251', '252', '253', '254', '255', '256', '257', '258', '259', '26', '260', '261', '262', '263', '264', '265', '266', '267', '268', '269', '27', '270', '271', '272', '273', '274', '275', '276', '277', '278', '279', '28', '280', '281', '282', '283', '284', '285', '286', '287', '288', '289', '29', '290', '291', '292', '293', '294', '295', '296', '297', '298', '299', '3', '30', '300', '301', '302', '303', '304', '305', '306', '307', '308', '309', '31', '310', '311', '312', '313', '314', '315', '316', '317', '318', '319', '32', '320', '321', '322', '323', '324', '325', '326', '327', '328', '329', '33', '330', '331', '332', '333', '334', '335', '336', '337', '338', '339', '34', '340', '341', '342', '343', '344', '345', '346', '347', '348', '349', '35', '350', '351', '352', '353', '354', '355', '356', '357', '358', '359', '36', '360', '361', '362', '363', '364', '365', '366', '367', '368', '369', '37', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '38', '380', '381', '382', '383', '384', '385', '386', '387', '388', '389', '39', '390', '391', '392', '393', '394', '395', '396', '397', '398', '399', '4', '40', '400', '41', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '6', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '7', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '8', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '9', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99'}}, ...
             'scoutfunc',  3, ...  % PCA        or % 1, ... % Mean
             'scouttime',  1, ...  % Before     or % 2, ...  % After
             'freqbands',  {sBandLabels{iBand}, [num2str(freqL.(sBandLabels{iBand})) ',' num2str(freqH.(sBandLabels{iBand}))], 'mean'}, ...
             'isorth',     1, ...
             'outputmode', 1);  % Save individual results (one file per input file)

         % Add a name tage
         sConnBandAEC.(sBandLabels{iBand}) = bst_process('CallProcess', 'process_add_tag', ...
             sConnBandAEC.(sBandLabels{iBand}), [], ...
             'tag',    sBandLabels{iBand}, ...
             'output', 2);  % Add to file name

         % Save Vertex x Vertex Conn Matrix: Each Freq Band
         connMatrix = file_fullpath(sConnBandAEC.(sBandLabels{iBand}).FileName);
         outpath = fullfile(hcp_dir, 'brainstormResults', ...
             'vertexAECConnectivity', SubjectNames{iSubj});
         if ~isfolder(outpath)
             mkdir(outpath)
             copyfile(connMatrix, fullfile(outpath, strcat(SubjectNames{iSubj}, ...
                 '_aecConn_', sBandLabels{iBand}, '.mat')))
         else
             copyfile(connMatrix, fullfile(outpath, strcat(SubjectNames{iSubj}, ...
                 '_aecConn_', sBandLabels{iBand}, '.mat')))
         end

        % Process: Phase locking value NxN
        sConnBandPLV.(sBandLabels{iBand}) = bst_process('CallProcess', 'process_plv1n', ...
        	sSrcRest, [], ...
            'timewindow', [sTime.Time(1) sTime.Time(end)], ...
            'scouts',     {'schaefer400.4k.label', {'0', '1', '10', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '11', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '12', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '13', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '14', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '15', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '16', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '17', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '18', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '19', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '2', '20', '200', '201', '202', '203', '204', '205', '206', '207', '208', '209', '21', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219', '22', '220', '221', '222', '223', '224', '225', '226', '227', '228', '229', '23', '230', '231', '232', '233', '234', '235', '236', '237', '238', '239', '24', '240', '241', '242', '243', '244', '245', '246', '247', '248', '249', '25', '250', '251', '252', '253', '254', '255', '256', '257', '258', '259', '26', '260', '261', '262', '263', '264', '265', '266', '267', '268', '269', '27', '270', '271', '272', '273', '274', '275', '276', '277', '278', '279', '28', '280', '281', '282', '283', '284', '285', '286', '287', '288', '289', '29', '290', '291', '292', '293', '294', '295', '296', '297', '298', '299', '3', '30', '300', '301', '302', '303', '304', '305', '306', '307', '308', '309', '31', '310', '311', '312', '313', '314', '315', '316', '317', '318', '319', '32', '320', '321', '322', '323', '324', '325', '326', '327', '328', '329', '33', '330', '331', '332', '333', '334', '335', '336', '337', '338', '339', '34', '340', '341', '342', '343', '344', '345', '346', '347', '348', '349', '35', '350', '351', '352', '353', '354', '355', '356', '357', '358', '359', '36', '360', '361', '362', '363', '364', '365', '366', '367', '368', '369', '37', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '38', '380', '381', '382', '383', '384', '385', '386', '387', '388', '389', '39', '390', '391', '392', '393', '394', '395', '396', '397', '398', '399', '4', '40', '400', '41', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '6', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '7', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '8', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '9', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99'}}, ...
            'scoutfunc',  3, ...  % PCA
            'scouttime',  1, ...  % Before     or % 2, ...  % After
            'freqbands',  {sBandLabels{iBand}, [num2str(freqL.(sBandLabels{iBand})) ',' num2str(freqH.(sBandLabels{iBand}))], 'mean'}, ...
            'mirror',     0, ...
            'keeptime',   0, ...
            'plvmeasure', 2, ...  % Magnitude
            'outputmode', 1);  % Save individual results (one file per input file)

        % Add a name tage
        sConnBandPLV.(sBandLabels{iBand}) = bst_process('CallProcess', 'process_add_tag', ...
            sConnBandPLV.(sBandLabels{iBand}), [], ...
            'tag',    sBandLabels{iBand}, ...
            'output', 2);  % Add to file name

        % Save Vertex x Vertex Conn Matrix: Each Freq Band
        connMatrix = file_fullpath(sConnBandPLV.(sBandLabels{iBand}).FileName);
        outpath = fullfile(hcp_dir, 'brainstormResults', ...
            'vertexPLVConnectivity', SubjectNames{iSubj});
        if ~isfolder(outpath)
            mkdir(outpath)
            copyfile(connMatrix, fullfile(outpath, strcat(SubjectNames{iSubj}, ...
                '_plvConn_', sBandLabels{iBand}, '.mat')))
        else
            copyfile(connMatrix, fullfile(outpath, strcat(SubjectNames{iSubj}, ...
                '_plvConn_', sBandLabels{iBand}, '.mat')))
        end


     end

    % Save report
    ReportFile = bst_report('Save', []);
    if ~isempty(reports_dir) && ~isempty(ReportFile)
        bst_report('Export', ReportFile, bst_fullfile(reports_dir, ...
            ['report_' ProtocolName '_' SubjectNames{iSubj} '.html']));
    end

%     clearvars -except SubjectNames BadChannels iSubj freqL freqH ProtocolName ...
%         reports_dir hcp_dir
    toc

end
