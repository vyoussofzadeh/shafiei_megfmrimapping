% Script generated by Brainstorm (19-Aug-2019)

% Input files
sFiles = {...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_02.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_03.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_04.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_05.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_06.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_07.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_08.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_09.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_10.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_11.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_12.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_13.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_14.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_15.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_16.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_17.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_18.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_19.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_20.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_21.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_22.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_23.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_24.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_25.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_26.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_27.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_28.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_29.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_30.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_31.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_32.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_33.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_34.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_35.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_36.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_37.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_38.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_39.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_40.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_41.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_42.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_43.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_44.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_45.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_46.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_47.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_48.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_49.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_50.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_51.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_52.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_53.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_54.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_55.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_56.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_57.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_58.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_59.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_60.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_61.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_62.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_63.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_64.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_65.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_66.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_67.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_68.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_69.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_70.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_71.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_72.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_73.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_74.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_75.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_76.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_77.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_78.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_79.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_80.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_81.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_82.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_83.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_84.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_85.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_86.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_87.mat', ...
    'corcoran_margaret/DFN_corcoran_margaret_IC_data/data_DFN_corcoran_margaret_IC_88.mat'};


% Start a new report
bst_report('Start', sFiles);

path_tools = '/usr/local/MATLAB_Tools';
set_ft
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions'));
addpath(genpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/run'));

% Process: Compute sources BF-FT
sFiles = bst_process('CallProcess', 'vy_ft_sourceanalysis_dics', sFiles, [], ...
     'Comment',     '', ...
    'method',      'dics', ...  % dSPM
    'mne',        struct(...
    'NoiseCov',      [], ...
    'InverseMethod', 'dics', ...
    'ChannelTypes',  {{}}, ...
    'SNR',           3, ...
    'diagnoise',     0, ...
    'SourceOrient',  {{'fixed'}}, ...
    'loose',         0.2, ...
    'depth',         1, ...
    'weightexp',     0.5, ...
    'weightlimit',   10, ...
    'regnoise',      1, ...
    'magreg',        0.1, ...
    'gradreg',       0.1, ...
    'eegreg',        0.1, ...
    'ecogreg',       0.1, ...
    'seegreg',       0.1, ...
    'fMRI',          [], ...
    'fMRIthresh',    [], ...
    'fMRIoff',       0.1, ...
    'pca',           1), ...
    'sensortypes', 'MEG, MEG MAG, MEG GRAD', ...
    'output',      1);  % Kernel only: shared


% Process: FieldTrip: ft_sourceanalysis
sFiles = bst_process('CallProcess', 'process_inverse_dics', sFiles, [], ...
    'method',     {'dics', {'LCMV beamformer', 'SAM beamformer', 'DICS beamformer', 'MNE', 'sLORETA', 'eLORETA', 'MUSIC', 'PCC', 'Residual variance'; 'lcmv', 'sam', 'dics', 'mne', 'sloreta', 'eloreta', 'music', 'pcc', 'rv'}}, ...
    'sensortype', {'MEG', {'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'; 'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'}});

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

