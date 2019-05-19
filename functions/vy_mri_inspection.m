function vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir, saveflag)

%%
figure;
ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_headshape(headshape);
ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
view ([0 90])
if isempty(saveflag)~=2
    savepath = fullfile(outputmridir,'headshape');
    hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
end
%%
% ft_determine_coordsys(individual_headmodel, 'interactive', 'no')
% hold on; % add the subsequent objects to the same figure
% ft_plot_headshape(headshape);
% view ([-10 40 10]);
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'headshape2');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end

%% headmodel - inspection
% figure;
% ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% hold on;
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% view ([-10 40 10]);
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'grid');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end

%%
% ft_determine_coordsys(mri_realigned, 'interactive', 'no')
% hold on;
% ft_plot_headshape(headshape);
% view ([50 80 10])
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% if isempty(saveflag)~=2
%     savepath = fullfile(outputmridir,'gridmri');
%     hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% end
% disp('---------------')
% disp(['figures were saved at,',outputmridir])

%%
% figure
% ft_plot_vol(individual_headmodel, 'unit', 'mm');  %this is the brain shaped head model volume
% ft_plot_sens(t_data.all.grad, 'unit', 'mm', 'coilsize', 10);  %this is the sensor locations  
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
% ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
