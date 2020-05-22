clear all; clc; close all;
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% =============
% Simulate eq0
% ============

% load and simulate
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
eq0 = read_eq(shot, time_ms/1000, cake_dir);

snow0 = analyzeSnowflake(eq0);
xp0 = [snow0.rx snow0.zx];
sfp = snow0.snowPlus;
sfm = ~sfp;
eq0 = designeq_ml(xp0,shot,time_ms);
sim0 = heatsim_ml(eq0,shot,time_ms);

% =========================
% Find new x-pts & simulate
% =========================
N = 10;

eqs{1}  = eq0;
sims{1} = sim0;
xps{1}  = xp0;


for k = 1:N
  fprintf(['\nIteration ' num2str(k) ' of ' num2str(N)])
  
  if sfm
    xps{k+1}  = estimate_xpts_sfm(eqs{k}, sims{k}, shot, time_ms, 1);
    eqs{k+1}  = designeq_ml(xps{k+1}, shot, time_ms, eqs{k});
    
  elseif sfp && ~constrain_sp
    xps{k+1} = estimate_xpts_sfp(eqs{k}, sims{k});
    eqs{k+1}  = designeq_ml(xps{k+1}, shot, time_ms, eqs{k});
    
  elseif sfp && constrain_sp
    sps = estimate_strike_pts(eqs{k}, sims{k}); 
    opts.constrain_sp = 1;
    opts.sp = sps;
    eqs{k+1} = designeq_ml(xp0, shot, time_ms, eqs{k}, opts);
    snf = analyzeSnowflake(eqs{k+1});
    xps{k+1} = [snf.rx snf.zx];                
  end
  
  sims{k+1} = heatsim_fit(eqs{k+1}, shot, time_ms, .006, .0025, .18, .18, 1);
%   sims{k+1} = heatsim_fit(eqs{k+1},shot,time_ms,1); 
end

eqs = {eqs{1}, eqs{end}};

save('xps','xps')
save('eqs','eqs')
save('sims','sims')
































