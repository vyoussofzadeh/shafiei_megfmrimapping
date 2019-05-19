function cln_data = vy_ica_cleaning_light(cfg_main, f_data)

n = cfg_main.n; % ICs

satis = 0;
disp('ica cleaning ...');
if exist(cfg_main.savepath, 'file') == 2
    load(cfg_main.savepath)
else
    
    comp = vy_ica(f_data,cfg_main.lay, n);
    title(savepath)
    cfg = [];
    cfg.updatesens = 'no';
    bic = input('Select bad ICs for slected data:');
    cfg.component = comp.label(bic);
    cln_data = ft_rejectcomponent(cfg, comp, f_data);
    close all
    if cfg_main.saveflag == 1
        save(cfg_main.savepath, 'cln_data', '-v7.3');
    end
end