clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked')

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/data/MEG/Clinical/MEG';
%- Output dir
outdir = '/data/MEG/Clinical';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

switch task
    case 1
        % - Auditory definition naming
        tag = 'DFN'; tag1 = 'dfn';
        Evnt_IDs = 1; % questions
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
end

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: 2016')
disp('5: 2015');
disp('6: 2014');
disp('7: 2013')
disp('8: 2012');
disp('9: 2011');
% disp('10: Older');
year = input('Year data were acquired: ');

clear ytag;
switch year
    case 1
        stag = {'19'};ytag = {'19'};
    case 2
        stag = {'18'};ytag = {'18'};
    case 3
        stag = {'17'};ytag = {'17'};
    case 4
        stag = {'16'};ytag = {'16'};
    case {5,6}
        stag = {'up'};ytag = {'14_5'};
%     case 6
%         stag = {'up'};ytag = {'14'};
    case 7
        stag = {'13'};ytag = {'13'};
    case 8
        stag = {'12'};ytag = {'12'};
    case 9
        stag = {'11'};ytag = {'11'};
end
disp('============');

%% Listing data
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
for j=1:numel(stag)
    stag1 = stag{1,j};
    d = rdir([indir,['/**/',stag1,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
    %     d = rdir([indir,['/**/',ytag1,'*/','sss','/*',tag,'*/*raw_sss.fif']]);
    if exist('tag1','var')
        d1 = rdir([indir,['/**/',stag1,'*/','sss','/*',tag1,'*/*raw_tsss.fif']]); d=[d;d1];
    end
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj_all{i} = [num2str(i), ': ', datafile{i}(Index(5)+1:Index(6)-1)];
    end
    datafile1 = vertcat(datafile1,datafile);
end
datafile1 = datafile1';
disp(datafile1)
disp('============');

%%
disp('1: choose specific subject');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        disp(datafile1)
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
        disp(subj_all')
end

%%
for i = 1:size(datafile2,1)
    
    datafile = datafile2{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Date:',Date])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    %     if year>=4
    %         yttag = 'older';
    %     else
    yttag = ytag{1};
    %     end
    outd.sub = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
%     pause(2)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    Run_preprocess
    
    %%
    if flag.notch ==1
        cfg = [];
        cfg.savefile = [];
        cfg.saveflag = 2;
        cfg.foilim = [2 40];
        cfg.plotflag  = 1;
        cfg.tapsmofrq = 5;
        cfg.taper     = 'hanning';
        [freq,ff,psd] = vy_fft(cfg, cln_data);
        grid on
        grid minor
        title(['Before band-stop filtering-',subj]);
        
        %- finding notch freq
        idx = find(freq.freq ==30); TF = islocalmax(psd); TF(1:idx-5) = 0;
        hold on
        plot(ff(TF),psd(TF),'r*')
               
        idx2 = find(TF == 1);
        if length(idx2)==1
            fsb = round(ff(idx2));
        else
            [val, idx] = max(psd(TF)); 
            fff = ff(TF); 
%             fsb = input(['Enter the sop-band frequency for ', subj,'?']);
            fsb = round(fff(idx));
        end
        disp([num2str(fsb),'Hz freq was selected for notch filteting']);
        
        %%
        cfg = [];
        cfg.bsfilter = 'yes';
        %     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
        cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
        %     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
        cln_data = ft_preprocessing(cfg, cln_data);
        %     cfg.bsfreq = [2 12]; % or whatever you deem appropriate
        
        cfg = [];
        cfg.savefile = [];
        cfg.saveflag = 2;
        cfg.foilim = [2 40];
        cfg.plotflag  = 1;
        cfg.tapsmofrq = 5;
        cfg.taper     = 'hanning';
        vy_fft(cfg, cln_data);
        grid on
        grid minor
        title(['After band-stop filtering-',subj]);
    end
    
    %% FFT & TFR                                                          
    if flag.freq == 1
        stag = 'tsk_baseline'; datain = cln_data;
        Run_freq
        disp(['time_of_interest:',num2str(time_of_interest),'sec']);
        disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
        L = 0.3;
    end
    
   %% Timelock  
    if flag.time == 1

        %--setting baseline interval
        toi(1,:) = [-0.3,0];
        
        %-- setting the post-stim interval
        disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
        %     disp('Yes: 1, No: 2');
        %     tfa = input('Is it OK to proceed?');
        %     disp('the following time was selected');
        if (time_of_interest-L) > 1.5 && (time_of_interest+L) < 2.3
            toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
        else
            switch task
                case 1
                    %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
                    toi(2,:) = [0.7,1.6]; % Best of DN
                case 2
%                     toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
                    toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
            end
        end
        disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
        disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
        Run_time
    end
    
    %% Grand Mean
    if flag.gave == 1
        Run_grandmean
    end
    %% Source analysis
    if flag.sourceanalysis == 1
        
        outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
        if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
        %     close all
        %% Processing       
        disp('1: Surface-based')
        disp('2: Volumetric');
        %     analysis = input('Eneter the analysis: ');
        analysis = 2;
        
        switch analysis
            case 1
                mtag = 'source_surf';
                disp('1: SPM source analysis (surface + BF)');
                disp('2: Export ft to Brainstorm');
                disp('3: Source-based bf');
                method = input('Method: ');
            case 2
                % end
                disp('1: LCMV')
                disp('2: Network+Connectvity');
                disp('3: DICS Source');
                %             method = input('Method: ');
                method = 3;
                %-
                %             disp('1: Low-res grid')
                %             disp('2: High-res grid')
                %             meshgrid = input('Mesh grid: ');
                meshgrid = 1;
        end
        
        switch analysis
            case 1
                % Surface-based analysis
                Run_surfacebased
            case 2
                % Volumetric-based analysis
                anatomy_check_flag = 2;
                Run_volumetric
        end
        disp('============');
    end
    %%
    pause(1)
    close all
    disp([datafile,' ,was completed'])
    disp(['output as,', outd.sub])
    
end
