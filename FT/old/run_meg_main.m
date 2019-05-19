clear; clc, close('all'); warning off

%% initial settings
restoredefaultpath
cd '\\utdrive.uthsc.edu\babajanilab\VNS\Scripts\MEG\ft';
cd_org = cd;
addpath(genpath(fullfile(cd_org,'functions')));
vy_init(cd_org)

load temp_grid_8mm % from, vy_warping()
hcp_path = 'F:\My Matlab\MEG\HCP\megconnectome-3.0';
template_mri = ft_read_mri(fullfile(hcp_path,'template','T1.nii')); %

DestDirectory = 'H:\VNS'; % saving directory

ft_old = 'F:\My Matlab\My codes\My GitHub\fieldtrip';

ft_path = 'F:\My Matlab\My codes\My GitHub\fieldtrip_packs\fieldtrip-20180809';
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal\ROI_MNI_V4.nii'));

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
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
end
datafolder = datafolder';
% disp(datafolder)

clear b
for i=1:length(datafolder)
    b{i} = num2str(i);
end
disp(table([b',datafolder]))


%%
sub = input('enter subject:');

%% details
% details_run_meg

%% step 1
run_meg_1st_SingleSubject
