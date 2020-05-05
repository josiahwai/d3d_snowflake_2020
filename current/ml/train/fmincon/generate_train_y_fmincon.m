% -------
% TESTING
% -------
clear all; clc; close all;
shot = 155353;
runbatch = 1;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));



% find times for which there are cake eqs and ir data
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
d = dir(cake_dir);
t_cake = [];
for k = 3:length(d)
  t_cake = [t_cake; str2num(d(k).name(end-3:end))];
end

qperp_dir  = [root 'inputs/qperp/']; 
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
t_ir = t;

[dt,k] = min(abs(t_cake - t_ir));  
times = t_cake(k(dt<20));  % within 20 ms

t =[    4391
        4358
        3528
        4291
        4325
        4491
        4524
        4424
        4590];
times(ismember(times,t)) = [];

% ----------------------
% GENERATE TRAIN Y
% ----------------------
for iTime = 1:length(times)
  try
  time_ms = times(iTime);  
  
  % load eq
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  cake_eq = read_eq(shot, time_ms/1000, cake_dir);
  cake_snow = analyzeSnowflake(cake_eq);
  xp0 = [cake_snow.rx cake_snow.zx];
    
  % -----------------------------------
  % GENERATE AND RUN BATCH JOBS
  % -----------------------------------    
  
  if runbatch
    
    % set up the batch jobs folders
    % .............................
    jobdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/'];       
    batchscript = [root 'ml/train/fmincon/sfd_fmincon_batch.sbatch'];    
    if exist(jobdir,'dir'), rmdir(jobdir,'s'); end
    mkdir(jobdir)            
    
    % copy x-pt initial condition to jobdir
    args = [shot time_ms xp0];
    
    save([jobdir '/args.mat'], 'args');
    
    % copy scripts to jobdir
    jobscript = [root 'ml/train/fmincon/sfd_fmincon_batch.m'];
    copyfile(jobscript, jobdir);
    
    % cd and submit batch job
    cd(jobdir)
    system(['sbatch ' batchscript]);
    cd(jobdir)
    
  end
  catch
  end
  cd([root 'ml/train'])
end








