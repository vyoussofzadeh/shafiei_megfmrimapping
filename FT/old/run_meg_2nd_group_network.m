clear; clc, close('all')

%% initial settings
vy_init

%% source averaging - whole brain
disp('1: CRM');
disp('2: VG-Auditory')
disp('3: VG-Printed')
task = input('Eneter the task number: ');
switch task
    case 1
        tsk = 'CRM'; % Continuous recognition memory
    case 2
        tsk = 'VGA'; % VerbGen-Auditory
    case 3
        tsk = 'VGP'; % VerbGen-Printed
end

disp('1: Source_lcmv');
disp('2: Source_lcmv-stats');
disp('3: Source_dics');
disp('4: Network_evc');
disp('5: SPM-COH(loreta)');
disp('6: all methods, lcmv, dics, coh, ');
method = input('Eneter the method: ');

% DestDirectory = fullfile('H:\VNS\Processed\MEG\ft3',tsk);
DestDirectory = fullfile('H:\VNS\Processed\MEG\ft',tsk);
d = dir(DestDirectory); files = {d.name}'; files_sel = files(3:end,1);

atlaspath = 'F:\My Matlab\My codes\My GitHub\fieldtrip_041718\fieldtrip-master\template\atlas';
atlas = ft_read_atlas(fullfile(atlaspath,'aal\ROI_MNI_V4.nii'));

switch method
    case 1
        mtd = 'lcmv';
        vy_group_source_lcmv_analysis
    case 2
        mtd = 'lcmv-stas';
        vy_group_source_lcmv_stats_analysis
    case 3
        mtd = 'dics';
        vy_group_source_dics_analysis
    case 4
        mtd = 'net';
        vy_group_network_analysis1
    case 5
        mtd = 'spm_coh';
        vy_group_source_spm_noica
    case 6
        mtd = 'all_sources';
%         vy_group_source_all
        vy_group_source_all_publication
%         vy_group_source_all_coh_fmri
end


