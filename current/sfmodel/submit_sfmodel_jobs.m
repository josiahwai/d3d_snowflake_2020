clear all; clc; close all;

% ========
% SETTINGS
% ========
shot = 155328;
runbatch = 1;
sim_sfm = 0;      % simulate times where IR predicts snowflake minus
sim_sfp = 1;      % simulate times where IR predicts snowflake plus
constrain_sp = 0; % for sfp only, constrain via strike pt instead of x-pt

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

% ir times
qperp_dir  = [root 'inputs/qperp/']; 
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t

t_ir = [];
if sim_sfm, t_ir = [t_ir t_3pks(1:20)]; end
if sim_sfp, t_ir = [t_ir t_2pks(1:20)]; end


[dt,k] = min(abs(t_cake - t_ir));  
t_sim = t_cake(k(dt<20));  % within 20 ms


% ============
% SUBMIT JOBS
% ============
if runbatch
for iTime = 1:length(t_sim)
  
  time_ms = t_sim(iTime);      
  
  batchscript_fn = [root 'sfmodel/sfmodel.sbatch'];   
  if sim_sfm  
    jobdir = [root 'sfmodel/jobs/' num2str(shot) '_sfm/' num2str(time_ms) '/'];
  elseif sim_sfp && constrain_sp
    jobdir = [root 'sfmodel/jobs/' num2str(shot) '_sfp_constrain_sp/' num2str(time_ms) '/'];    
  elseif sim_sfp && ~constrain_sp
    jobdir = [root 'sfmodel/jobs/' num2str(shot) '_sfp/' num2str(time_ms) '/'];
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

  cd([root 'sfmodel'])
end
end





