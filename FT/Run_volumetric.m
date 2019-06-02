%% Volumetric-based analysis
mridir = fullfile(indir,subj,'brainstorm_db/anat');
d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
clear fid
if ~isempty(d)
    sMRI1 = d.name;
    load(sMRI1);
    fid.SCS = SCS;
    fid.NCS = NCS;
    mripfile = fullfile(mridir,'T1.nii');
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    
    cfg = [];
    cfg.megdata = t_data.app;
    cfg.mripfile = mripfile;
    cfg.hsfile = datafile; % headshape;
    cfg.fid = fid;
    cfg.outputmridir = outputmridir;
    cfg.subj = subj;
    cfg.plotflag = 2;
    [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg);
    %     [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm, individual_seg] = vy_mri_neuromag4(cfg);
    %         [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag3(cfg);
    
end

%% Choosing mesh
switch meshgrid
    case 1
        meshtag = 'lowres';
        load temp_grid % low-res
        template_grid = ft_convert_units(template_grid, 'mm');
        individual_grid = individual_grid_10mm;
    case 2
        meshtag = 'highres';
        load temp_grid_8mm % high-res
        individual_grid = individual_grid_8mm;
end

%% Anatomoy check!
saveflag = 2;
if anatomy_check_flag == 1
    vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
end
%     close all

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

%%
switch method
    case 1
        %%
        vy_source_lcmv_light
        %                             vy_source_lcmv
    case 2
        %%
        mtag = 'conn'; outd.vol = fullfile(outd.sub,mtag);
        cfg = [];
        cfg.p.ft_old = allpath.ft_old;
        cfg.p.ft_path = allpath.ft_path;
        cfg.p.cd_org = cd_org;
        cfg.p.hcp_path = allpath.hcp_path;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        vy_network_light1(cfg,t_data) % conn-network analysis
        
    case 3
        %%
        mtag = 'dics'; outd.vol = fullfile(outd.sub,mtag);
        cfg = [];
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.sens = sens;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        vy_source_dics(cfg,ep_data);
        
        %     case 4
        %         %%
        %         outputdir = fullfile(outdir,'ft_process',yttag, subj, tag);
        %         outputdir1 = fullfile(outputdir, 'spm_source');
        %         if exist(outputdir1, 'file') == 0
        %             mkdir(outputdir1);   % create a directory
        %         end
        %         cfg = [];
        %         cfg.toilim = [-0.4 2];
        %         eint_data = ft_redefinetrial(cfg, cln_data);
        %
        %         if exist(mripfile, 'file') == 2
        %             cd(outputdir1);
        %
        %             cfg = [];
        %             cfg.p.spm = spm_path;
        %             cfg.p.hcp_path = hcp_path;
        %             cfg.p.ft_path = ft_path;
        %             cfg.p.spmbf = spmbf_path;
        %             cfg.p.cd_org = cd_org;
        %             cfg.datafile = datafile;
        %             cfg.eint_data = eint_data;
        %             cfg.mripfile = mripfile;
        %             cfg.subj = subj;
        %             vy_forward_spm_meg(cfg);
        %
        %             restoredefaultpath
        %             addpath(genpath(ft_path));
        %             addpath(genpath(hcp_path));
        %             addpath(genpath([cd_org,'/functions']));
        %             addpath(genpath([cd_org,'/Data_file']));
        %             cd(outputdir)
end
