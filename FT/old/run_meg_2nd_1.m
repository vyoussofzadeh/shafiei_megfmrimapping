clear; clc, close('all')

%% initial settings
vy_init

%% source averaging - whole brain
disp('1: CRM');
disp('2: VG-Auditory')
task = input('Eneter the task number: ');
switch task
    case 1
        tsk = 'CRM'; % Continuous recognition memory
    case 2
        tsk = 'VGA'; % VerbGen-Auditory
    case 3
        tsk = 'VGP'; % VerbGen-Printed
end

DestDirectory = fullfile('H:\VNS\Processed\MEG\ft',tsk);
d = dir(DestDirectory); files = {d.name}'; files_sel = files(3:end,1);

%% Group network analysis
vy_group_network


