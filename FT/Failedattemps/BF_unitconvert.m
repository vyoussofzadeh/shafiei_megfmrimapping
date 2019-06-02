clear
sourcemodel = ft_read_headshape('tess_cortex_mid_low.mat');

sMri = '/data/MEG/Clinical/MEG/bednar_peggy/brainstorm_db/anat/bednar_peggy/brainstormsubject.mat';

save('sourcemodel.mat','sourcemodel');

% ROI coordinates
for i = 1:length(sourcemodel.pos)
%     A = At.Scouts(i).Seed;
    newPosScs = sourcemodel.pos(i,:);
    newPosVox = round(cs_convert(sMri, 'scs', 'voxel', newPosScs));
    newPosMNI(i,:) = round((1e3.*cs_convert(sMri, 'voxel', 'mni', newPosVox)),1);
    roi{i}= At.Scouts(i).Region;
    roi_l{i}= At.Scouts(i).Label;
end


%%
transform1 = inv(neuromag_channels.TransfMeg{strcmp('neuromag_device=>neuromag_head',neuromag_channels.TransfMegLabels)});
transform2 = inv(neuromag_channels.TransfMeg{strcmp('neuromag_head=>scs',neuromag_channels.TransfMegLabels)});
transform3 = inv(neuromag_channels.TransfMeg{strcmp('refine registration: head points',neuromag_channels.TransfMegLabels)});


%%
channels2 = Channels;
channels2.Channel=channels2.Channel(1:306);

clear a b
a = cln_data.label;
for i=1:length(channels2.Channel)
    b{i,:} = channels2.Channel(i).Name;
end
[sharedvals,idx] = intersect(b,a,'stable');
channels2.Channel=channels2.Channel(idx);
channels2.Projector = [];


