
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
d = rdir(fullfile(bsdatadir,subj,['/@*',tag,'*'],'/channel_vectorview306*.mat'));

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

%%
save(fullfile(bsdatadir,'BS_channels'),'chan_updated');
save(fullfile(bsdatadir,'IC_data'),'cln_data');


