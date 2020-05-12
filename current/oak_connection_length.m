clear all; clc; close all
root = '/u/jwai/d3d_snowflake_2020/current/';
shot = 165288;
time_ms = 4200;

% =============
% Simulate eq0
% =============

% load and simulate
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
eq0 = read_eq(shot, time_ms/1000, cake_dir);
struct_to_ws(eq0.gdata);
load('d3d_obj_mks_struct_129129.mat');
struct_to_ws(tok_data_struct);

plot_eq(eq0)
contour(rg,zg,psizr,psibry + [-0.1:.01:0])
r0 = 2.285;
z0 = 0;

[L, r1, z1] = calcConnectionLength(r0, z0, psizr, rg, zg, ...
    bzero, rzero, limdata, 1,1)
  
scatter(r0,z0,'k','filled')
scatter(r1,z1,'k','filled')








