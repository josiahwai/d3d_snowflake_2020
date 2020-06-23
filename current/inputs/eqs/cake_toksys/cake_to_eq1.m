clear all; clc; close all; 
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

% =========
% SETTINGS
% =========

shotlist = 155336:2:155344;
saveit = 1;

% =============================
% Convert cake to toksys format
% =============================

topdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/cake/';
save_topdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/cake_toksys/';

for shot = shotlist
  
  cake_dir = [topdir num2str(shot) '/'];
  d = dir(cake_dir);
  d(1:2) = [];
  
  for k = 1:length(d)
    try
    
      time_ms = str2double( d(k).name(end-3:end));
      cake_eq = read_eq(shot, time_ms/1000, cake_dir);
      cake_snow = analyzeSnowflake(cake_eq);
      xp = [cake_snow.rx cake_snow.zx];
      
      eq = designeq_ml(xp,shot,time_ms);
      eq = rmfield(eq, {'r', 'b', 'p'});
      
      if saveit
        savedir = [save_topdir num2str(shot) '/'];
        if ~exist(savedir, 'dir')
          mkdir(savedir)
        end
        fn = ['eq' num2str(shot) '_' num2str(time_ms)];
        save([savedir fn], 'eq')
      end
      
    catch 
    end
  end
end
  
