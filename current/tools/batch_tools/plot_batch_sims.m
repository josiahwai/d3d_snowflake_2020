clear all; clc; close all;
% ========
% SETTINGS
% ========
shotdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/155355_sfm';


% ===================
% LOAD AND PLOT SIMS
% ===================
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

shotdir_info = dir(shotdir);

nJobs = length(shotdir_info) - 2;

for k = 1:nJobs
  try
    time_ms = str2num(shotdir_info(k+2).name)
    
    jobdir = [shotdir '/' num2str(time_ms) '/'];
    
    load([jobdir 'eqs.mat'])
    load([jobdir 'sims.mat'])
    load([jobdir 'xps.mat'])
    
    analyze_sfmodel_sim
  catch
  end  
end

