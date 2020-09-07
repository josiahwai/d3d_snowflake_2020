% load('eq0')
% shot = 155340;
% time_ms = 3694;
% 
% xp = [1.1677 1.12178 -1.0973 -1.2559];
% % xp = [1.166 1.23 -1.10 -1.25];
% 
% % eq  = designeq_dum(xp, shot, time_ms, eq0);
% 
% 
% [rxP,rxS,zxP,zxS] = unpack(xp);
% spec.locks.rx = [rxP rxS];
% spec.locks.zx = [zxP zxS];
% load('d3d_obj_mks_struct_6565.mat')
% config = tok_data_struct;
% init = eq0;
% 
% eq = gsdesign(spec, init, config);

eqdir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_snowtypes/149743/gEQDSK';

eq = read_eq(149743, 3.9, eqdir);

plot_eq(eq)




