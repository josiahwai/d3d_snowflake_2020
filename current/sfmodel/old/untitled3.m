clear all; clc; close all
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root))

shot = 155354;

load('/u/jwai/d3d_snowflake_2020/current/inputs/qperp/qperp_155354.mat')

for k = 1:length(t_2pks)
  time_ms = floor(t_2pks(k));
  
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  eq0 = read_eq(shot, time_ms/1000, cake_dir);
  
  figure(9)
  clf
  plot_eq(eq0)
  title(num2str(time_ms))
  axis([1 1.5 -1.4 -.9])
end


