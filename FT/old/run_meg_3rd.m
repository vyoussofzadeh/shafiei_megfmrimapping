clear; clc, close('all')

%% initial settings
vy_init

%% source averaging - whole brain
disp('1: CRM');
disp('2: VG-Auditory')
disp('3: VG-Printed')
task = input('Enter the task number: ');
switch task
    case 1
        tsk = 'CRM'; % Continuous recognition memory
    case 2
        tsk = 'VGA'; % VerbGen-Auditory
    case 3
        tsk = 'VGP'; % VerbGen-Printed
end

DestDirectory = 'H:\VNS\Processed\MEG\ft'; % saving directory
outputdir = fullfile(DestDirectory,tsk, 'group');
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end

disp('1: network vs fmri');
disp('2: source (spm_noica) vs fmri');
disp('3: source (dics_noica) vs fmri');
disp('4: source (coh_noica) vs fmri');
method = input('Enter the method: ');
switch method
    case 1
        load(['.\Data_file\fmri_', tsk,'.mat']);
        d = rdir([outputdir,'\net\indiv_n_wb\rhighres_*.nii']);
        clear datafolder
        for i=1:length(d)
            datafolder{i} = d(i).name;
        end
        datafolder = datafolder';
        disp('MEG:')
        disp(datafolder);
        savelabel = 'net';
    case 2
        load(['.\Data_file\fmri_', tsk,'.mat']);
        load(['.\Data_file\meg_spm_', tsk,'.mat']);
    case 3
        load(['.\Data_file\fmri_', tsk,'.mat']);
        load(['.\Data_file\meg_dics_12_23Hz_', tsk,'.mat']);
        disp(data_dis)
    case 4
        load(['.\Data_file\fmri_', tsk,'.mat']);
        load(['.\Data_file\meg_coh_', tsk,'.mat']);
        disp(data_dis)
end
disp('fMRI')
disp(niifile);

%%
cd(outputdir)

%% fmri - nii
niifile_fmri = niifile;

clear vol_fmri
for i=1:size(niifile_fmri,1)
    [fmri_dir,b] = fileparts(niifile_fmri(i,:));
    niifile1 = fullfile(fmri_dir,[b,'.nii']);
    s_vol_fmri = ft_read_mri(niifile1);
    vol_fmri(i,:,:,:) = s_vol_fmri.anatomy;
end

vol_name_fmri = fullfile(fmri_dir,'group',[tsk,'.nii']);
switch task
    case 1
        projthresh_fmri = 0.4;
    case 2
        projthresh_fmri = 0.2;
    case 3
        projthresh_fmri = 0.3;
end

s_vol_fmri.anatomy = vy_normalize(vol_fmri);
s_vol_fmri = vy_vol_thresh(s_vol_fmri,projthresh_fmri); % abs
s_vol_fmri1 = vy_vol_thresh_posneg(s_vol_fmri,.001); % keep positive effects only!

Opt = [];
Opt.savenii = 1;
Opt.savefig = 1;
Opt.savename = fullfile(outputdir,['fmri_',tsk]);
vy_surfce_vis2(s_vol_fmri1,vol_name_fmri, Opt);

%-
[~, par_fmri_group, ~] = vy_parcellate(s_vol_fmri1, atlas,'anatomy');
save('par_fmri_group','par_fmri_group')
% load par_fmri_group

%% meg - nii
niifile_meg = datafolder;
meg_dir = fileparts(niifile_meg{1});

clear vol_meg
for i=1:numel(niifile_meg)
    [meg_dir,b] = fileparts(niifile_meg{i});
    niifile1 = fullfile(meg_dir,[b,'.nii']);
    s_vol_meg = ft_read_mri(niifile1);
    s_vol_meg.anatomy(isnan(s_vol_meg.anatomy(:)))=0;
    vol_meg(i,:,:,:) = s_vol_meg.anatomy;
end

savedir = outputdir;
vol_name_meg = fullfile(savedir,[tsk,'.nii']);

projthresh_meg = 0.6;
s_vol_fmri1 = vy_vol_thresh_posneg(s_vol_fmri,.001); % keep positive effects only!
s_vol_meg.anatomy = squeeze(mean(vol_meg,1));

s_vol_fmri1.anatomy(isnan(s_vol_fmri1.anatomy(:)))=0;

a = (s_vol_fmri1.anatomy); a = a./max(a(:));
b = (s_vol_meg.anatomy); b = b./max(b(:));
% save('s_vol_meg','s_vol_meg')
load s_vol_meg

s_vol_meg = vy_vol_thresh(s_vol_meg,projthresh_meg); % abs
s_vol_meg_group = vy_vol_thresh_posneg(s_vol_meg,.001); % keep positive effects only!

Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = fullfile(savedir,['meg_',tsk]);
vy_surfce_vis2(s_vol_meg_group,vol_name_meg, Opt);

%-
[~, par_meg_group, ~] = vy_parcellate(s_vol_meg_group, atlas,'anatomy');
save(['par_meg_group_',tsk],'par_meg_group')
% load (['par_meg_group_',tsk])

%% plot parcellation (meg)
addpath(spm_path);
addpath(connpath);

sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = 'anatomy';
cfg1.downsample   = 1;
% cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, par_meg_group, ft_read_mri(sMRI, 'format', 'nifti_spm'));

savedir = fullfile(meg_dir,tsk,'group');

Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = fullfile(savedir,'Parc_meg');
vy_surfce_vis2(sourceint_pow,fullfile(savedir,'Parc_meg.nii'), Opt);

%% Masking combined meg-fmri
%-inspired by, https://www.nitrc.org/forum/forum.php?thread_id=7218&forum_id=1144
savemask = fullfile(savedir,['meg_net_fmri_',tsk,'.nii']);
spm_imcalc({vol_name_meg, vol_name_fmri},savemask, 'i1.*(i1>0)-i2.*(i2>0)');

Opt = [];
Opt.savenii = 2;
Opt.savefig = 1;
Opt.savename = fullfile(savedir,['meg_net_fmri_',tsk]);
vy_surfce_vis2([],savemask, Opt);

%% Individuals (fmri)
savedir = fullfile(fmri_dir, 'indv');
if exist(savedir, 'file') == 0, mkdir(savedir); end

%-fmri
clear vol_name_fmri_indv
projthresh_fmri = 0.3;
for i = 1:size(niifile_fmri,1)
    
    s_vol_fmri.anatomy = squeeze(vol_fmri(i,:,:,:));
    s_vol_fmri = vy_vol_thresh(s_vol_fmri,projthresh_fmri); % abs
    s_vol_fmri1 = vy_vol_thresh_posneg(s_vol_fmri,.001,1); % keep positive effects only!
    
    nii = niifile_fmri(i,:);
    Index = strfind(nii, '_');
    switch task
        case 1
            subj_fmri{i} = nii(Index(end-1)+1:Index(end)-1);
        case 2
            subj_fmri{i} = nii(Index(end-2)+1:Index(end-1)-1);
        case 3
            subj_fmri{i} = nii(Index(end-2)+1:Index(end-1)-1);
    end
    vol_name_fmri_indv{i} = fullfile(fmri_dir,'indv',[tsk, '_', subj_fmri{i},'.nii']);
    
%     %-Surf-vis
    opt.run = [];
    opt.tsk = tsk;
    opt.subj = subj_fmri{i};
    opt.savedir = savedir;
    opt.savenii = 1;
    vy_surfce_vis(s_vol_fmri1,vol_name_fmri_indv{i}, opt);
    
    %     [~, par_fmri, ~] = vy_parcellate(s_vol_fmri1, atlas,'anatomy');
    %     par_fmri_indv(i,:) = par_fmri.anatomy;
end
% save('par_fmri','par_fmri_indv','par_fmri','subj_fmri','vol_name_fmri_indv')
load('par_fmri');

%-parcellation check
name = 'fmri';
par_fmri.anatomy = mean(par_fmri_indv,1)';
vy_parcellate_plot(par_fmri, coor, name);

%%
%-meg
clear vol_name_meg_ind par_meg_indv
projthresh_meg = 0.4;
savedir = fullfile(meg_dir, tsk, 'group','indiv');
if exist(savedir, 'file') == 0, mkdir(savedir); end
for i=1:size(niifile_meg,1)
    
    s_vol_meg.anatomy = squeeze(vol_meg(i,:,:,:));
    %     s_vol_meg = vy_vol_thresh(s_vol_meg,projthresh_meg); % abs
    %     s_vol_meg1 = vy_vol_thresh_posneg(s_vol_meg,.001,1); % keep positive effects only!
    
    nii = niifile_meg(i,:);
    clear Index
    Index = strfind(nii{1}, '_');
    switch method
        case 1
            subj_meg{i} = nii{1}(Index(end-1)+1:Index(end)-1);
            run = nii{1}(Index(end)+1:end-4);
        case 3
            subj_meg{i} = nii(Index(end-1)+1:Index(end)-1);
            run = nii(Index(end)+1);
        case 4
            subj_meg{i} = nii(Index(4)+1:Index(5)-1);
            run = nii(Index(5)+1:Index(6)-1);
    end
    vol_name_meg_ind{i} = fullfile(savedir,[tsk, '_', subj_meg{i},'_',run,'.nii']);
    
    a = (s_vol_meg_group.anatomy);a(isnan(a(:)))=0; a = a./max(a(:));
    b = (s_vol_meg.anatomy); b(isnan(b(:)))=0; b = b./max(b(:));
    s_vol_meg1 = vy_vol_thresh(s_vol_meg1,projthresh_meg); % abs
    s_vol_meg2 = vy_vol_thresh_posneg(s_vol_meg1,.001, 1); % keep positive effects only!
    %
    %     Opt = [];
    %     Opt.savenii = 1; Opt.savefig = [];
    %     Opt.savename = fullfile(savedir,['meg_',tsk]);
    %     vy_surfce_vis2(s_vol_meg2,vol_name_meg_ind{i}, Opt);
    
    %     %-Surf-vis
    %     Opt = [];
    %     opt.run = run;
    %     opt.tsk = tsk;
    %     opt.subj = subj_meg{i};
    %     opt.savedir = savedir;
    %     opt.savenii = 1;
    %     opt.savefig = 1;
    %     vy_surfce_vis(s_vol_meg2,vol_name_meg_ind{i}, opt);
    %
    %     %- save nii only
    %     cfg = [];
    %     cfg.filetype  = 'nifti';
    %     cfg.datatype   = 'uint8'; %'float';
    %     cfg.parameter = 'anatomy';
    %     cfg.filename  = vol_name_meg_ind{i};
    %     ft_volumewrite(cfg, s_vol_meg2);
    
    [~, par_meg, ~] = vy_parcellate(s_vol_meg1, atlas,'anatomy');
    par_meg_indv(i,:) = par_meg.anatomy;
    
    %     [~, par_meg, ~] = vy_parcellate(s_vol_meg2, atlas,'anatomy');
    %     par_meg_indv(i,:) = par_meg.anatomy;
end
save('par_meg_indv','par_meg_indv','par_meg','vol_name_meg_ind')
load('par_meg_indv');

%-parcellation check
name = 'net';
par_meg.anatomy = mean(par_meg_indv,1)';
vy_parcellate_plot(par_meg, coor, name);

%%
%- Overlay meg-fmri
projthresh_meg = 0.4;
k=1;
clear par_meg_indv
Index = strfind(meg_dir, '\');
savedir = fullfile(meg_dir,tsk,'group','meg_fmri');
if exist(savedir, 'file') == 0, mkdir(savedir); end
for i=1:length(subj_meg)
    for j=1:length(subj_fmri)
        if strncmp(subj_fmri{j},subj_meg{i},5)==1
            
            disp([subj_fmri{j},'-vs-',subj_meg{i}])
            
            f = vol_name_fmri_indv{j};
            Index1 = strfind(f, '\');  Index2 = strfind(f, '.nii');
            savemask = fullfile(savedir,['comb_',f(Index1(end)+1:Index2-1),'.nii']);
            
            %             s_vol_fmri = vy_vol_thresh_posneg(s_vol_fmri,.001, 1); % keep positive effects only!
            s_vol_fmri = ft_read_mri(vol_name_fmri_indv{j});
            s_vol_fmri.anatomy(isnan(s_vol_fmri.anatomy(:)))=0;
            
            
            s_vol_meg = s_vol_fmri;
            s_vol_meg.anatomy = squeeze(vol_meg(i,:,:,:));
            
            a = (s_vol_fmri1.anatomy); a = a./max(a(:));
            b = (s_vol_meg.anatomy); b = b./max(b(:));
            s_vol_meg = vy_vol_thresh(s_vol_meg,projthresh_meg); % abs
            
            %             vol_name_meg{i} = fullfile(savedir,[tsk, '_', subj_meg{i},'_',run,'.nii']);
            
%             Opt = [];
%             Opt.savenii = 1; Opt.savefig = [];
%             Opt.savename = fullfile(savedir,['meg_',tsk]);
%             vy_surfce_vis2(s_vol_meg,vol_name_meg_ind{i}, Opt);
            
            nii = vol_name_meg_ind{i};
            clear Index
            Index = strfind(nii, '_');
            run = nii(Index(end)+1:end-4);
                        
            opt = [];
            opt.run = run;
            opt.tsk = tsk;
            opt.subj = subj_meg{i};
            opt.savedir = savedir;
            opt.savenii = 1;
            opt.savefig = 1;
            vy_surfce_vis(s_vol_meg,vol_name_meg_ind{i}, opt);
            
            spm_imcalc({vol_name_meg_ind{i}, vol_name_fmri_indv{j}},savemask, 'i1.*(i1>0)-i2.*(i2>0)');
            
            %-Surf-vis
            opt.run = run;
            opt.tsk = ['comb-',tsk];
            opt.subj = subj_fmri{j};
            opt.savedir = savedir;
            opt.savenii = 2;
            vy_surfce_vis([],savemask, opt);
            masked{k} = subj_fmri{j}; k=k+1;
                       
            %-parcellation (meg)
            [~, par_meg1, ~]  = vy_parcellate(s_vol_meg, atlas,'anatomy');
            par_meg_comb(i,:) = par_meg1.anatomy;
            %             vy_parcellate_plot(par_meg1, coor, name);            
        end
    end
end
save('par_meg_comb','par_meg_comb','subj_meg')
% load('par_meg_comb');

%% Plot laterality (template)
idx = [12:2:20,80:2:88];
idxx = [idx, idx-1];
par_meg1 = par_meg;
par_meg1.anatomy = zeros(size(par_meg.anatomy,1),1);
% for i=1:length(idxx)
%     par_meg.anatomy(idxx(i)) = rand(1);
% end
for i=1:length(idxx)
    par_meg1.anatomy(idxx(i)) = 3.*rand(1);
end
% vy_parcellate_plot(par_meg, coor, 'lat_template');

sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = 'anatomy';
cfg1.downsample   = 1;
cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, par_meg1, ft_read_mri(sMRI, 'format', 'nifti_spm'));

savedir = fullfile(meg_dir,tsk,'group');
Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = fullfile(savedir,'lat_template');
vy_surfce_vis2(sourceint_pow,fullfile(savedir,'lat_template.nii'), Opt);
% manual colormap
% cmap = vega20c(2*96);
% state.colormap=cmap;

%% Correlation between meg-net and fMRI
fmri = par_fmri.anatomy./max(par_fmri.anatomy(:));
fmri_sd = std(par_fmri_indv)';
meg = par_meg.anatomy./max(par_meg.anatomy(:));
meg_sd = (std(par_meg_indv)./sqrt(size(par_meg_indv,1)))';
% save('par_meg_indv1','meg','meg_sd');
load par_meg_indv1

% fmri = par_fmri.anatomy;
% meg = par_meg.anatomy;

figure,bar([fmri,meg])
[~,p] = corr(fmri,meg);
legend([{'fmri'};{'meg-network'}])
L = length(meg);
set(gca,'Xtick', 1:L,'XtickLabel',par_meg.label);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','arial');
grid
box off
set(gca,'color','none');
set(gcf, 'Position', [500   100   2000   500]);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','Arial');

%%
idx = [12:2:20,80:2:88];
idx2 = [idx, idx-1];

% par_meg.anatomy = meg2;
% sourceint_pow = vy_parcellate_plot(par_meg, coor, 'net');
bar_input=[fmri(idx2),meg(idx2)]';
errorbar_input=[fmri_sd(idx2),meg_sd(idx2)]';
label = {'Frontal-Inf-Oper-R' ;
    'Frontal-Inf-Tri-R'  ;
    'Frontal-Inf-Orb-R'  ;
    'Rolandic-Oper-R'    ;
    'Supp-Motor-Area-R'  ;
    'Heschl-R'           ;
    'Temporal-Sup-R'     ;
    'Temporal-Pole-Sup-R';
    'Temporal-Mid-R'     ;
    'Temporal-Pole-Mid-R';
    'Frontal-Inf-Oper-L' ;
    'Frontal-Inf-Tri-L'  ;
    'Frontal-Inf-Orb-L'  ;
    'Rolandic-Oper-L'    ;
    'Supp-Motor-Area-L'  ;
    'Heschl-L'           ;
    'Temporal-Sup-L'     ;
    'Temporal-Pole-Sup-L';
    'Temporal-Mid-L'     ;
    'Temporal-Pole-Mid-L'};
errorbar_groups(bar_input,errorbar_input, 'bar_names',label,'bar_width',0.75,'errorbar_width',0.5);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','arial');
% grid
box off
set(gca,'color','none');
set(gcf, 'Position', [500   500   1000   400]);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gca,'FontName','Arial');
legend([{'fmri'};{'meg-network'}]);
r = corr2(bar_input(1,:),bar_input(2,:));

%%
%-meg
D = [];
D.label = label';
D.eigenvector_cent = 4.*bar_input(2,:)';
[ROI, ROI_sel] = vy_ROI_report(D,0.4,coor(idx2,:), 'eigenvector_cent');

%-fmri
D = [];
D.label = label';
D.tstats = 4.*bar_input(1,:)';
[ROI, ROI_sel] = vy_ROI_report(D,0.6,coor(idx2,:), 'tstats');


%==============================================
%- plot parcellation (meg)
meg3 = meg;
meg3((meg < 0.8.*max(meg(:)))) = nan;
par_meg.anatomy = meg3;
sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = 'anatomy';
cfg1.downsample   = 1;
% cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, par_meg, ft_read_mri(sMRI, 'format', 'nifti_spm'));

savedir = fullfile(meg_dir,tsk,'group');
Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = fullfile(savedir,'Parc_meg');
vy_surfce_vis2(sourceint_pow,fullfile(savedir,'Parc_meg.nii'), Opt);
view([-114,28])
view([114,28])

%- plot parcellation (fmri)
fmri2 = fmri;
fmri2((fmri2 < 0.7.*max(fmri2(:)))) = nan;
par_meg.anatomy = fmri2;
sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = 'anatomy';
cfg1.downsample   = 1;
% cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, par_meg, ft_read_mri(sMRI, 'format', 'nifti_spm'));

savedir = fullfile(meg_dir,tsk,'group');
Opt = [];
Opt.savenii = 1; Opt.savefig = [];
Opt.savename = fullfile(savedir,'Parc_fmri');
vy_surfce_vis2(sourceint_pow,fullfile(savedir,'Parc_fmri.nii'), Opt);

%% Laterality index
close all
load laterality_CRM
figure(1),
subplot 121
hbar = bar(par_crm);
view([90 -90])
L = length(par_crm);
for i=1:length(idx)
    S{i} = num2str(idx(i));
end
% set(gca,'Xtick',idx,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
title('Word-recognition','fontsize',16)
% legend([{'fMRI-tValue'};{'MEG'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);

DataArray = par_crm;
Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
figure(2)
subplot 121
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(par_crm,'Label',{'fMRI','meg'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
title('Word-recognition','fontsize',16)
disp(['fmri: ',num2str(mean(par_crm(:,1))), '+-', num2str(std(par_crm(:,1)))]);
disp(['meg: ',num2str(mean(par_crm(:,2))), '+-', num2str(std(par_crm(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end


for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end

%
load laterality_VGA
figure(1)
subplot 122
hbar = bar(par_vga);
view([90 -90])
L = length(par_vga);
for i=1:length(idx)
    S{i} = num2str(idx(i));
end
% set(gca,'Xtick',idx,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
title('Verb-Generation','fontsize',16)
% legend([{'fmri-tValue'};{'meg-hubs'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
disp(['fmri: ',num2str(mean(par_vga(:,1))), '+-', num2str(std(par_vga(:,1)))]);
disp(['meg: ',num2str(mean(par_vga(:,2))), '+-', num2str(std(par_vga(:,2)))]);
% colormap(Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);


for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end

figure(2)
subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(par_vga,'Label',{'fmri','meg'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
set(gca,'color','none');
set(gca,'FontName','HelveticaNeueLT Std Lt');
title('Verb-Generation','fontsize',16)
xlabel('Modality','FontSize', 16);
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end
