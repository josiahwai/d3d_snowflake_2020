% -------
% TESTING
% -------
clear all; clc; close all;
shot = 165288;
runbatch = 1;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

efit_dir = [root 'inputs/eqs/cake/' num2str(shot)];
d = dir(efit_dir);
times = [];
for k = 3:length(d)
  times = [times; str2num(d(k).name(end-3:end))];
end

% ----------------------
% GENERATE TRAIN Y
% ----------------------
for iTime = 1:length(times)
  try
  time_ms = times(iTime);  
  
  % load eq
  efit_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  efit_eq = read_eq(shot, time_ms/1000, efit_dir);
  efit_snow = analyzeSnowflake(efit_eq);
  xp = [efit_snow.rx efit_snow.zx];
    
  % -----------------------------------
  % GENERATE AND RUN BATCH JOBS
  % -----------------------------------
  
  % 2 different initial conditions centered around xp
  dxp = .01 * [0 -1 1 0; 0 1 -1 0];
  
  if runbatch
    
    % set up the batch jobs folders
    % .............................
    job_topdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/'];
    output_dir = [root 'ml/train/job_outputs/'];
    batchscript = [root 'ml/train/sfd_fmincon_batch.sbatch'];
    
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    if exist(job_topdir,'dir'), rmdir(job_topdir,'s'); end
    mkdir(job_topdir)
    
    
    
    for iJob = 1:1
      % copy x-pt initial condition to jobdir
      xp0 = xp + dxp(iJob,:);
      args = [iJob shot time_ms xp0];
      
      jobdir = [job_topdir num2str(iJob)];
      mkdir(jobdir);
      save([jobdir '/args.mat'], 'args');
      
      % copy scripts to jobdir
      jobscript = [root 'ml/train/sfd_fmincon_batch.m'];
      copyfile(jobscript, jobdir);
      
      % cd and submit batch job
      cd(jobdir)
      system(['sbatch ' batchscript]);
      cd(job_topdir)
    end
    
    cd([root 'ml/train'])
  end
  catch
  end
end








