clear; clc; close all;
plotit = 1;
shot = 165288;
time_ms = 4200;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath('/u/jwai/d3d_snowflake_2020/current')); 
load('/u/jwai/d3d_snowflake_2020/current/inputs/tok_data/d3d_obj_mks_struct_6565.mat')
efitdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/efit01/165288';

eq = read_eq(shot,time_ms,efitdir);
struct_to_ws(eq.gdata);

limdata = tok_data_struct.limdata;
spRZ = isoflux_spFinder(psizr, psibry, rg, zg, limdata, 1:length(limdata));


plot_eq(eq)
axis([1.0 1.5 -1.4 -0.9])
scatter(spRZ(:,1),spRZ(:,2),'g','filled')


[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 500, true);
psilim = bicubicHermite(rg,zg,psizr,rlim,zlim);

