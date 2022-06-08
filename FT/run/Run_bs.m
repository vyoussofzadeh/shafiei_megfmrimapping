
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
d = rdir(fullfile(bsdatadir,subj,['/@*',tag,'*'],'/channel_vectorview306*.mat'));
if isempty(d)
    d = rdir(fullfile(bsdatadir,subj,['/*',tag,'*'],'/channel_vectorview306*.mat'));
    if exist('tag1','var') && isempty(d)
        d = rdir(fullfile(bsdatadir,subj,['/*',tag1,'*'],'/channel_vectorview306*.mat'));
    end
end

%%
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    %     databf{i} = pathstr;
    datafilebf{i} = d(i).name;
end
disp(datafilebf')
channel = load(datafilebf{1});
disp('=============')
disp(datafilebf{1})
disp('was loaded')

%%
% 5: Editing channel file: run below script (while BS is open)


%%
bssavedir = fullfile(indir,subj,'brainstorm_db/');
switch task
    case {1,2}
        chan_updated = channel;
        chan_updated.Channel=chan_updated.Channel(1:306);
        clear a b
        a = cln_data.label;
        for i=1:length(chan_updated.Channel)
            b{i,:} = chan_updated.Channel(i).Name;
        end
        [sharedvals,idx] = intersect(b,a,'stable');
        chan_updated.Channel=chan_updated.Channel(idx);
        chan_updated.Projector = [];
        disp('=============')
        disp('channel modification was completed');
        save(fullfile(bssavedir,[tag,'_',subj,'_BS_channels']),'chan_updated');
        save(fullfile(bssavedir,[tag,'_',subj,'_IC_data']),'cln_data','-v7.3');
    case 3
        clear a b chan_updated
        chan_updated = channel;
        chan_updated.Channel=chan_updated.Channel(1:306);
        %--- left
        cln_data = cln.left;
        a = cln_data.label;
        for i=1:length(chan_updated.Channel)
            b{i,:} = chan_updated.Channel(i).Name;
        end
        [sharedvals,idx] = intersect(b,a,'stable');
        chan_updated.Channel=chan_updated.Channel(idx);
        chan_updated.Projector = [];
        save(fullfile(bssavedir,[tag,'_',subj,'_left_BS_channels']),'chan_updated');
        %         save(fullfile(bssavedir,[tag,'_',subj,'_left_IC_data']),'cln_data','-v7.3');
        
        clear a b chan_updated
        chan_updated = channel;
        chan_updated.Channel=chan_updated.Channel(1:306);
        %-- right
        cln_data = cln.right;
        a = cln_data.label;
        for i=1:length(chan_updated.Channel)
            b{i,:} = chan_updated.Channel(i).Name;
        end
        [sharedvals,idx] = intersect(b,a,'stable');
        chan_updated.Channel=chan_updated.Channel(idx);
        chan_updated.Projector = [];
        %         save(fullfile(bssavedir,[tag,'_',subj,'_right_IC_data']),'cln_data','-v7.3');
        save(fullfile(bssavedir,[tag,'_',subj,'_right_BS_channels']),'chan_updated');
end

%% Converting to fif file
% isepch = ft_datatype(cln_data, 'raw');
% hdr = ft_read_header(datafile);
% test_var = cln_data;
% test_var.hdr = hdr;
%
% if exist('Run') ==2
%     fieldtrip2fiff(fullfile(bssavedir,[tag,'_',subj,'_IC_data.fif']), test_var)
% else
%     fieldtrip2fiff(fullfile(bssavedir,[tag,'_',subj,'_',Run,'_IC_data.fif']), test_var)
% end

%%
restoredefaultpath
addpath((allpath.ft_path));
ft_defaults
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));
