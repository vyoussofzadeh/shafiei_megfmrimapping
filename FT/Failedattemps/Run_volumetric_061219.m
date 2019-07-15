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
    cfg.atlas = atlas;
    [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm, mri_realigned_ctf] = vy_mri_neuromag2(cfg);
end

%% Choosing mesh
switch meshgrid
    case 1
        meshtag = 'lowres';
        %         load('standard_sourcemodel3d10mm');
        load temp_grid % low-res
        template_grid = ft_convert_units(template_grid, 'mm');
        individual_grid = individual_grid_10mm;
    case 2
        meshtag = 'highres';
        %         load('standard_sourcemodel3d8mm');
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
        %% LCMV
%         mtag = 'lcmv';
        mtag = 'lcmv_bl';
%         mtag = 'lcmv_stat';
        outd.vol = fullfile(outd.sub,mtag);        
        cfg = [];
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.mtag = mtag;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.toi       = toi;
        cfg.fb   = [12,30]; % Hz
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        cfg.mtd = mtag;
        vy_source_lcmv(cfg,cln_data);

%     case 2
%         %% lcmv-stat, beta
%         mtag = 'lcmv_stat'; outd.vol = fullfile(outd.sub,mtag);
%         cfg = [];
%         cfg.allpath = allpath;
%         cfg.grid = individual_grid;
%         cfg.headmodel = individual_headmodel;
%         cfg.subj = subj;
%         cfg.sens = sens;
%         cfg.atlas = allpath.atlas_path;
%         cfg.outputdir = outd.vol;
%         cfg.toi       = toi;
%         cfg.template_grid = template_grid;
%         cfg.template_mri = template_mri;
%         cfg.mtd = mtag;
%         vy_source_lcmv_beta(cfg,cln_data);
%         %         vy_source_lcmv_light
%         %                             vy_source_lcmv
    case 3
        %%
        mtag = 'conn'; 
        mtag = 'conn_bl';
        outd.vol = fullfile(outd.sub,mtag);
        cfg = [];
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.mtag = mtag;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        vy_network_light1(cfg,t_data) % conn-network analysis
%     case 4
%         %%
%         mtag = 'conn_bl'; outd.vol = fullfile(outd.sub,mtag);
%         cfg = [];
%         cfg.allpath = allpath;
%         cfg.grid = individual_grid;
%         cfg.headmodel = individual_headmodel;
%         cfg.subj = subj;
%         cfg.sens = sens;
%         cfg.mtag = mtag;
%         cfg.fb   = [12,30]; % Hz
%         cfg.toi  = toi;
%         cfg.atlas = allpath.atlas_path;
%         cfg.outputdir = outd.vol;
%         cfg.template_grid = template_grid;
%         cfg.template_mri = template_mri;
%         vy_network_light1(cfg,cln_data) % conn-network analysis
        
    case 5
        %%
        mtag = 'dics';
        mtag = 'dics_stat';
        outd.vol = fullfile(outd.sub,mtag);
        cfg = [];
        cfg.grid = individual_grid;
        cfg.allpath = allpath;
        cfg.headmodel = individual_headmodel;
        cfg.sens = sens;
        cfg.subj = subj;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        vy_source_dics(cfg,ep_data);
        
%     case 6
%         %% DICS stats
%         mtag = 'dics_stat'; outd.vol = fullfile(outd.sub,mtag);
%         cfg = [];
%         cfg.grid = individual_grid;
%         cfg.allpath = allpath;
%         cfg.headmodel = individual_headmodel;
%         cfg.sens = sens;
%         cfg.subj = subj;
%         cfg.outputdir = outd.vol;
%         cfg.template_grid = template_grid;
%         cfg.template_mri = template_mri;
%         %         cfg.effect  = effects;
%         vy_source_dics_stats(cfg,ep_data);
        
end
