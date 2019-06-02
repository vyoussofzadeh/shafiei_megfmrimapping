clear; clc, close('all'); warning off

brainstorm
addpath(genpath('./functions'));

%%
% datafile = '/data/MEG/Vahab/test_data/raghavan_manoj2/brainstorm_db/anat/raghavan_manoj2';
% datafile = '/data/MEG/Vahab/Scripts/Vahab/Test_MEG2/bednar_peggy/brainstorm_db/anat/bednar_peggy';
datafile = '/data/MEG/Clinical/MEG/dougherty_danielle/brainstorm_db/anat/dougherty_danielle';

disp(datafile)
d = rdir(fullfile(datafile,'subjectimage*.mat'));
if ~isempty(d)
    sMRI = d.name;
    cd(datafile)
    cd ..
    OutputFile = fullfile(pwd,'T1.nii');
    out_mri_nii(sMRI, OutputFile);
end
