% clear all; clc; close all;
warning('on', 'all');
% ========
% SETTINGS
% ========

% shotdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm_efit/155354_sfm/';
shotdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/';

cd(shotdir)

% ===================
% LOAD AND PLOT SIMS
% ===================
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

shotdir_info = dir(shotdir);

nJobs = length(shotdir_info) - 2;
% tgood = [];

for k = 1:nJobs
  try
    time_ms = str2num(shotdir_info(k+2).name);
    
    jobdir = [shotdir '/' num2str(time_ms) '/'];
    
%     % check if ran
%     outfile = dir([jobdir '*.out']);
%     fid = fopen(outfile.name);
%     nIter = 0;
%     tline = fgetl(fid);    
%     while ischar(tline) 
%       if contains(tline, 'Iteration')
%         nIter = nIter + 1;
%       end
%       tline = fgetl(fid);
%     end
%     fclose(fid);

    % plot stuff
    if ismember(time_ms, tgood)
      
%       tgood = [tgood; time_ms];      
      load([jobdir 'sims.mat'])
      load([jobdir 'xps.mat'])
      load([jobdir 'eqs.mat'])
      
      close all
      plotsim(sims{1})
      plotsim(sims{end})
      title(num2str(time_ms))
      
      figure(28)
      hold on
      psin = eqs{2}.psibar;
      plot(psin, eqs{2}.jpar, 'r', 'linewidth', 2)
      plot(psin, eqs{3}.jpar, 'b',  'linewidth', 2)
      xlim([0. 1.02])
      
      % analyze_sfmodel_sim
      
    else
      warning(['Did not run ' num2str(time_ms)])
    end    
    
  catch
  end  
end

