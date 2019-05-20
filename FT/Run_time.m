%% Epoching
%     switch task
%         case 1 %'DefNam'
%             toi = input('Eneter toi (e.g. [-0.3,0;1.5,2]): ');
%         case 2 %'PicNam'
%             %             toi = [-0.3,0;0.7,1];
%             toi = [-0.2,0;0.5,0.8];
%     end
%     toi = [-0.2,0;1,1.5];
% toi = [-0.3,0;0.7,1.5];
toi = [-0.3,0;0,1.5];
ep_data = vy_epoch(cln_data, toi);

%- Appending data
cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

%% Data-covarinace estimation
t_data = vy_timelock(ep_data);
% savepath = fullfile(outputdir,['tl_',subj,'_',run,'.mat']);
% save(savepath, 't_data', '-v7.3');
