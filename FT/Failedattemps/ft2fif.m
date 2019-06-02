
mne_path = '/data/MEG/Vahab/Github/MCW-MEGlab/tools/ft_packages/fieldtrip_master/external/mne';
addpath(mne_path);

fiff_write_evoked('test','cln_data')


%%
restoredefaultpath
ft19_path = fullfile(allpath.path_tools,'/ft_packages/fieldtrip_20190419');
addpath(ft19_path);
ft_defaults

%%
test_var = ep_data.all;
test_var = t_data.pst;
fieldtrip2fiff('test.fif', test_var)

hdr = ft_read_header(datafile);

test_var = cln_data;
test_var.hdr = hdr;
fieldtrip2fiff('test.fif', cln_data)

%%
test_var = ft_preprocessing([],cln_data);


%%
data   = ft_checkdata(test_var, 'datatype', {'raw', 'timelock'}, 'feedback', 'yes');

istlck = ft_datatype(data, 'timelock')
isepch = ft_datatype(data, 'raw')

%%
fieldtrip2fiff('test.fif', test_var)

%%
fieldtrip2fiff('test.fif', f_data)

%%
fiff_write_evoked('test.fif', test_var);

%%
% using older ver of ft for network analysis
restoredefaultpath
addpath(genpath(allpath.ft_old));
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));

%%
save('test.mat','f_data');



