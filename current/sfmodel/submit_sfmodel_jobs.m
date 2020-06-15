function submit_sfmodel_jobs(shot, sim_sfm, sim_sfp, constrain_sp)

% clear all; clc; close all;

% ========
% SETTINGS
% ========
% shot = 155330;
runbatch = 1;
% sim_sfm = 0;      % simulate times where IR predicts snowflake minus
% sim_sfp = 1;      % simulate times where IR predicts snowflake plus
% constrain_sp = 1; % for sfp only, constrain via strike pt instead of x-pt

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

% ===========
% FIND TIMES
% ===========

% Find times for the shot where: 
% cake eq exists, ir predicts wanted sfp/sfm
% ..........................................

% times with cake
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
d = dir(cake_dir);
t_cake = [];
for k = 3:length(d)
  t_cake = [t_cake; str2num(d(k).name(end-3:end))];
end

% times with efit
efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

d = dir(cake_dir);
t_efit = [];
for k = 3:length(d)
  t_efit = [t_efit; str2num(d(k).name(end-3:end))];
end


% ir times
qperp_dir  = [root 'inputs/qperp/']; 
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t

t_ir = [];
if sim_sfm, t_ir = [t_ir t_3pks]; end
if sim_sfp, t_ir = [t_ir t_2pks]; end

% simulate times with ir data, efit, and cake all within 10ms 
t_eq = intersect(t_cake, t_efit);

[dt,k] = min(abs(t_eq - t_ir)); 
t_sim = t_eq(k(dt<10));  % within 10 ms


% ============
% SUBMIT JOBS
% ============
if runbatch
for iTime = 1:length(t_sim)
  
  time_ms = t_sim(iTime);      
  
  batchscript_fn = [root 'sfmodel/sfmodel.sbatch'];   
  if sim_sfm  
    jobdir = [root 'sfmodel/jobs/sfm/' num2str(shot) '_sfm/' num2str(time_ms) '/'];
  elseif sim_sfp && constrain_sp
    jobdir = [root 'sfmodel/jobs/sfp_constrain_sp/' num2str(shot) '_sfp_sp/' num2str(time_ms) '/'];    
  elseif sim_sfp && ~constrain_sp
    jobdir = [root 'sfmodel/jobs/sfp/' num2str(shot) '_sfp/' num2str(time_ms) '/'];
  end
  
  
  if exist(jobdir,'dir'), rmdir(jobdir,'s'); end
  mkdir(jobdir)
  
  args = [shot time_ms constrain_sp];
  save([jobdir 'args.mat'], 'args');
  
  % copy scripts to jobdir
  jobscript = [root 'sfmodel/sfmodel.m'];
  copyfile(jobscript, jobdir);
  
  % cd and submit batch job
  cd(jobdir)
  system(['sbatch ' batchscript_fn]);
  cd(jobdir)

%   cd([root 'sfmodel'])
end
end





