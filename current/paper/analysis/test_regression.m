clear

load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp.mat')


struct_to_ws(sim);

% dist = sqrt((dxp(:,1) - dxp(:,3)).^2 + (dxp(:,2) - dxp(:,4)).^2);
% regression(dist, djmax, 1, 1)


regression(ddpsi, djmax, 1, 0)







