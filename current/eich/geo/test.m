clear; clc; close all

load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155328_sfm/3594/eqs.mat')
load('/u/jwai/d3d_snowflake_2020/current/inputs/qperp/qperp_155328.mat')
load('d3d_obj_mks_struct_6565.mat')

eq = eqs{end};
[~,k] = min(abs(t-3594));
qperp = qperp(k,:)';
 
ef = eich_fitter_geo(s, qperp, eq, tok_data_struct, 1)

ef = eich_fitter(s, qperp, eq, tok_data_struct, 1)


plot_eq(eq)





