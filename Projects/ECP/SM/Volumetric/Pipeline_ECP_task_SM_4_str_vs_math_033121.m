clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.trialinfo.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 0;
flag.freq = 0;     % TFR & FFT
flag.toi = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.speech = 0;             % speech analysis
flag.sourceanalysis = 1;     % source analysis
flag.anatomy = 1;            % grand average analysis
flag.anatomy_check = 0;
flag.meshgrid_sel = 1; % high-res = 2, low-res = 1;

%% Initial settings
% set(0,'DefaultFigureWindowStyle','normal');
% set(gcf,'units','points','position',[500,500,500,500]);
% set(0, 'DefaultFigureRenderer', 'painters');
% set(gcf, 'renderer', 'painters');
% set(0, 'DefaultFigureRenderer', 'opengl');

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/group/jbinder/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/data/MEG/Research/ECP/';

ECP_scriptdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';
% ECP_datadir = '/data/MEG/Research/ECP/';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
%- Story math
tag = [];
tag.task = 'SM';

%% Listing data
clc
% if exist(['datalog_',tag.task,'.mat'], 'file') == 2
%     load(['datalog_',tag.task,'.mat'])
% else
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
d = rdir([indir,['/**/tSSS/*',tag.task,'*_raw.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
end
datafile1 = vertcat(datafile1,datafile);
datafile1 = datafile1';

save(['datalog_',tag.task,'.mat'],'datafile1','subj_all')
% end
disp(datafile1)
disp('============');

%%
disp('1: choose specific subject');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        disp(subj_all')
        subsel1 = input('enter subject number?');
end
disp('============');

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%%
clear datafile2
switch subsel
    case 1
        datafile2{1} = datafile1{subsel1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    case 2
        datafile2 = datafile1;
end

%%
disp('1: Surface-based')
disp('2: Volumetric');
analysis = input('Eneter the analysis: ');
%     analysis = 2;
switch analysis
    case 1
        mtag = 'source_surf';
        disp('1: Surface-based SPM source analysis');
        disp('2: Export ft to Brainstorm');
        disp('3: Surface-based bf');
        method = input('Method: ');
    case 2
        % end
        disp('1: LCMV')
        disp('2: Network+Connectvity');
        disp('3: DICS Source');
        method = input('Method: ');
        %             method = 3;
        %-
        %             disp('1: Low-res grid')
        %             disp('2: High-res grid')
        %             meshgrid_sel = input('Mesh grid: ');
        meshgrid_sel = 1;
end

%%
disp('1: Str')
disp('2: Math');
disp('3: str_math');
dc = input('Select data condition:');

switch dc
    case 1
        tag.dcon = 'str';
    case 2
        tag.dcon = 'math';
    case 3
        tag.dcon = 'str_math';
end

%%
for dd = 1:size(datafile2,1)
    
    datafile = datafile2{dd}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(6)+1:Index(7)-1);
    Index = strfind(datafile, 'run');
    run  = datafile(Index+3);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Run:',run])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    outd.sub = fullfile(outdir,'MEG_work_ft',subj, tag.task, run);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    ic_selection = 1; % 1: manual, 2: automated
    Run_preprocess_SM_3
    
    %%
    idx_str = trlInfoColDescr1(:,2)==1;
    idx_math = trlInfoColDescr1(:,2)==2;
    
    %%
    %     %-Narration interval start, 8
    %     [C,IA,IB] = intersect(trlInfoColDescr1(:,8),trlInfoColDescr1(:,20));
    %
    %     %-Option interval start, 11
    %     [C,IA,IB] = intersect(trlInfoColDescr1(:,11),trlInfoColDescr1(:,20));
    
    %% Response time
    %     savefile = ['RT_',subj,'.mat'];
    %     if exist(savefile, 'file') == 2
    %         load(savefile);
    %     else
    %
    %         RT = (trlInfoColDescr1(:,17) - trlInfoColDescr1(:,15))./2e3;
    %         RT_str = RT(idx_str); idx_str2 = RT_str > 0; RT_str_uniq = unique(RT_str(idx_str2));
    %         RT_math = RT(idx_math); idx_math2 = RT_math > 0; RT_math_uniq = unique(RT_math(idx_math2));
    %
    %         m_RT_str_uniq = mean(RT_str_uniq);
    %         m_RT_math_uniq = mean(RT_math_uniq);
    %
    %         disp(['RT of str was: sec', num2str(m_RT_str_uniq)]);
    %         disp(['RT of math was: sec', num2str(m_RT_math_uniq)]);
    %
    %         save(savefile, 'RT', 'm_RT_str_uniq', 'm_RT_math_uniq');
    %     end
    
    %% Difficulty level
    %     savefile = ['DL_',subj,'.mat'];
    %     %     if exist(savefile, 'file') == 2
    %     %         load(savefile);
    %     %     else
    %     Run_triggers_SM
    %     %      load(['trial_info_',subj,'.mat']);
    %     idx_str = trlInfoColDescr(:,2)==1;
    %     idx_math = trlInfoColDescr(:,2)==2;
    %
    %     DL = (trlInfoColDescr(:,7));
    %     DL_str = DL(idx_str); idx_str2 = find(DL_str > 0); DL_str_uniq = unique(DL_str(idx_str2));
    %     DL_math = DL(idx_math); idx_math2 = find(DL_math > 0); DL_math_uniq = unique(DL_math(idx_math2));
    %
    %     DL_str_uniq
    %     DL_math_uniq
    %
    %     save(savefile, 'DL', 'DL_str_uniq', 'DL_math_uniq');
    %     end
    %% Accuracy
    %     savefile = ['Acc_',subj,'.mat'];
    %     if exist(savefile, 'file') == 2
    %         load(savefile);
    %     else
    %         idx_str = find(trlInfoColDescr1(:,2)==1);
    %         idx_math = find(trlInfoColDescr1(:,2)==2);
    %
    %         accuracy = length(find(trlInfoColDescr1(:,18)==1))/length(trlInfoColDescr1)*100;
    %         accuracy_str  = length(find(trlInfoColDescr1(idx_str,18)==1))/length(idx_str)*100;
    %         accuracy_math = length(find(trlInfoColDescr1(idx_math,18)==1))/length(idx_math)*100;
    %
    %         disp(['Accuracy was: %', num2str(accuracy)]);
    %         disp(['Accuracy of str was: %', num2str(accuracy_str)]);
    %         disp(['Accuracy of math was: %', num2str(accuracy_math)]);
    %
    %         save(savefile, 'accuracy', 'accuracy_str', 'accuracy_math');
    %     end
    %%
    %     length(find(trlInfoColDescr1(:,21) == 10)) % 10.Story Sentence
    %     length(find(trlInfoColDescr1(:,21) == 12)) % 12.Story Option Intro
    %     length(find(trlInfoColDescr1(:,21) == 13))  % 13.Math Option 1 Word
    %     length(find(trlInfoColDescr1(:,21) == 14))  % 14. Math Option OR Word
    %     length(find(trlInfoColDescr1(:,21) == 15))  % 15.Math Option 2 Word
    % %
    %     length(find(trlInfoColDescr1(:,21) == 20)) % 20.Math Narration number word
    %     length(find(trlInfoColDescr1(:,21) == 21)) % 21.Math Narration operand word
    %     length(find(trlInfoColDescr1(:,21) == 22)) % 22.Math Option Intro Wor
    %     length(find(trlInfoColDescr1(:,21) == 23)) % 23.Math Option 1 Word
    %     length(find(trlInfoColDescr1(:,21) == 24)) % 24. Math Option OR Word
    %     length(find(trlInfoColDescr1(:,21) == 25)) % 25.Math Option 2 Word
    %
    %'21 Event Type - 20.Math Narration number word,\n 21.Math Narration operand word,\n 22.Math Option Intro Word,\n 23.Math Option 1 Word\n  24. Math Option OR Word\n 25.Math Option 2 Word\n 10.Story Sentence,\n 12.Story Option Intro\n 13.Math Option 1 Word\n  14. Math Option OR Word\n 15.Math Option 2 Word'
    
    %     if length(cln_data.trialinfo) == length(trlInfoColDescr1)
    
    %         idx_str = find(trlInfoColDescr1(:,2)==1);
    %         idx_math = find(trlInfoColDescr1(:,2)==2);
    
    idx_str = find(cln_data.trialinfo(:,21)==10);
    idx_math = find(cln_data.trialinfo(:,21)==20);
    
    cfg = [];
    cfg.trials = idx_str;
    cln_data_str = ft_selectdata(cfg,cln_data);
    cfg.trials = idx_math;
    cln_data_math = ft_selectdata(cfg,cln_data);
    %     else
    %         error('something is wrong, preprocess data again!');
    %     end
    
    %% Notch filtering of 30Hz
    if flag.notch ==1
        cln_data = cln_data_str; Run_notch;cln_data_str = cln_data;
        cln_data = cln_data_math; Run_notch; cln_data_math = cln_data;
    end
    close all,
    
    %% FFT & TFR
    fmax = 40;
    if flag.freq == 1
        switch dc
            case 1
                datain = cln_data_str;
                tag.dcon = 'str';
                %                 savetag = fullfile(savepath,[stag,'_tfr_str_',subj,'.mat']);
                
                %                 time_of_interest_str = time_of_interest;
                %                 freq_of_interest_str = freq_of_interest;
                %                 disp(['time_of_interest:',num2str(time_of_interest_str),'sec']);
                %                 disp(['freq_of_interest:',num2str(freq_of_interest_str),'Hz']);
            case 2
                datain = cln_data_math;
                tag.dcon = 'math';
                %                 %                 savetag = fullfile(savepath,[stag,'_tfr_math_',subj,'.mat']);
                %                 time_of_interest_math = time_of_interest;
                %                 freq_of_interest_math = freq_of_interest;
                %                 disp(['time_of_interest:',num2str(time_of_interest_math),'sec']);
                %                 disp(['freq_of_interest:',num2str(freq_of_interest_math),'Hz']);
            case 3
                %                 datain1 = cln_data_str; datain2 = cln_data_math;
        end
        %         Run_freq
    end
    L = 0.3;
    
    if flag.sourceanalysis == 1
        
        %% Epoching
        %         switch dc
        %             case 1
        %                 toi_str = [-0.4,0;1.5,2.5];
        %                 %             toi_str = [-0.3,0;1.4,2];
        %                 %             toi_str = [-0.3,0;1,2];
        %                 %             toi_str = [-0.3,0;1,1.8];
        %                 ep_data_stor = vy_epoch(cln_data_str, toi_str);
        %                 stag = ['Str_',run];
        %
        %             case 2
        %                 toi_math = [-0.3,0;0.5,1.2];
        %                 ep_data_math = vy_epoch(cln_data_math, toi_math);
        %                 stag = ['Math_',run];
        %             case 3
        %                 %     toi_str_math = [0.4,1.5];
        %                 toi_str_math = [0,3];
        %                 %     toi_str_math = [0.5,1];
        %                 ep_data_str_math = vy_epoch_paired(cln_data_math, cln_data_str, toi_str_math);
        %         end
        %
        %%
        if flag.toi == 1
            %             toi(1,:) = [-0.3,0];
            %             toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
            %
            %             disp('===========================================');
            %             disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),';',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec was selected as contrast intervals']);
            %             %         disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as post-stim']);
            %             %         disp('===========================================');
            %             warning(['Maximum trial length:[', num2str(datain.time{1}(1)), ',', num2str(datain.time{1}(end)),']']);
            %             disp('OK to proceed: 1, No, another time interval: 2:');
            %             ask_time = input(':');
            %             if ask_time == 2
            %                 disp('Enter time interval in sec, eg, [-0.3,0; 0.7,1.2];');
            %                 toi = input(':');
            %             end
            
            %%
            switch dc
                case 1
                    datain = cln_data_str;
                    tag.dcon = 'str';
                 case 2
                    datain = cln_data_math;
                    tag.dcon = 'math';
            end
            
            for i=1:length(datain.trial)
                
                disp([num2str(i),'/', num2str(length(datain.trial))]);
                
                cfg         = [];
                cfg.trials = i;
                data_sel = ft_redefinetrial(cfg, datain);
                
                cfg         = [];
                cfg.length  = 1;
                cfg.overlap = 0;
                data_sel_seg = ft_redefinetrial(cfg, data_sel);
                
                cfg         = [];
                cfg.trials = 1:2;
%                 cfg.toilim = [-1,0];
                data_bsl = ft_redefinetrial(cfg, data_sel_seg);
                
                cfg         = [];
                cfg.trials = 3:length(data_sel_seg.trial);
                data_pst = ft_redefinetrial(cfg, data_sel_seg);
                
                if i==1
                    data_merged_bsl = data_bsl;
                    data_merged_pst = data_pst;
                else
                    cfg = [];
                    data_merged_bsl =  ft_appenddata(cfg, data_merged_bsl, data_bsl);
                    cfg = [];
                    data_merged_pst =  ft_appenddata(cfg, data_merged_pst, data_pst);
                end
                
            end
            %%
            ep_data = [];
            ep_data.bsl = data_merged_bsl;
            ep_data.pst = data_merged_pst;
            toi = [1,1;1,1];
            
            %%
            %%
            switch dc
                case 1
                    ep_data_stor = ep_data;
                    stag = ['Str_',run];
                    
                case 2
                    ep_data_math = ep_data;
                    stag = ['Math_',run];
                case 3
                    toi_str_math = [0,3];
                    ep_data_str_math = vy_epoch_paired(cln_data_math, cln_data_str, toi_str_math);
            end
            
            %% Appending data
            switch dc
                case 1
                    cfg = [];
                    ep_data_stor.app = ft_appenddata(cfg,ep_data_stor.bsl,ep_data_stor.pst);
                case 2
                    cfg = [];
                    ep_data_math.app = ft_appenddata(cfg,ep_data_math.bsl,ep_data_math.pst);
                case 3
                    a = length(ep_data_str_math.pst.label); b = length(ep_data_str_math.bsl.label);
                    [C,ia,ib] = intersect(ep_data_str_math.bsl.label,ep_data_str_math.pst.label, 'stable');
                    cfg = [];
                    if b > a
                        cfg.channel = ia;
                        ep_data_str_math.bsl = ft_selectdata(cfg, ep_data_str_math.bsl);
                    else
                        cfg.channel = ib;
                        ep_data_str_math.pst = ft_selectdata(cfg, ep_data_str_math.pst);
                    end
                    ep_data_str_math.app = ft_appenddata(cfg,ep_data_str_math.bsl,ep_data_str_math.pst);
            end
        end
        %% Timelock
        %     t_data_stor = vy_timelock(ep_data_stor);
        %     t_data_math = vy_timelock(ep_data_math);
        
        %% Grand Mean
        if flag.gave == 1
            datain = cln_data_str;  savefile = fullfile(savepath,['gmean_str_',subj,'.mat']); Run_grandmean_SM
            datain = cln_data_math; savefile = fullfile(savepath,['gmean_math_',subj,'.mat']); Run_grandmean_SM
        end
        
        %% Source analysis
        outputmridir = fullfile(outdir,'MEG_work_ft', subj,'anat'); % output dir
        if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
        
        %% Processing
        switch analysis
            case 1
                switch dc
                    case 1
                        datain = cln_data_str;
                    case 2
                        datain = cln_data_math;
                end
                bsdir = '/rcc/stor1/projects/ECP/MEG/MEG_work_BS';
                bsdatadir_help = '/rcc/stor1/projects/ECP/MEG/MEG_work_BS/BAK/data';
                % Surface-based analysis
                Run_surfacebased_ECP
            case 2
                % Volumetric-based analysis
                rl = 2;
                %             Run_volumetric
                freq_of_interest = 20; fmax=  40;
                switch dc
                    case 1
                        ep_data = ep_data_stor; flag.savetag = 'stor';
                    case 2
                        ep_data = ep_data_math; flag. savetag = 'math';
                    case 3
                        ep_data = ep_data_str_math; toi = toi_str_math; flag. savetag = 'stor_math';
                end
                Run_volumetric_ECP
        end
        disp('============');
        
        %
        h =  findobj('type','figure'); n = length(h);
        if n > 10
            close all
        end
        
        disp([datafile,' ,was completed'])
    end
end