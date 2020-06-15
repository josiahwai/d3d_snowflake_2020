clear all; clc; close all;
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% =============
% Simulate eq0
% ============

% load and simulate
efit_dir = '/p/omfit/users/jwai/projects/155354_efits/EFITtime/OUTPUTS/gEQDSK';

efit_eq1 = read_eq(shot, time_ms/1000, efit_dir);
efit_sim1 = heatsim_ml(efit_eq1.gdata,shot,time_ms);
efit_snow1 = analyzeSnowflake(efit_eq1);
efit_xp1 = [efit_snow1.rx efit_snow1.zx];

efit_xp2 = efit_xp1;
efit_eq2 = designeq_ml(efit_xp2,shot,time_ms);
efit_sim2 = heatsim_ml(efit_eq2,shot,time_ms);
efit_snow2 = analyzeSnowflake(efit_eq2);
  
  
% =========================
% Find new x-pts & simulate
% =========================
N = 15;

eqs = {efit_eq1, efit_eq2};
sims = {efit_sim1, efit_sim2};
xps  = {efit_xp1, efit_xp2};

sfp = efit_snow2.snowPlus;
sfm = ~sfp;



for k = 2:N
  fprintf(['\nIteration ' num2str(k) ' of ' num2str(N)])
  
  if sfm
    xps{k+1}  = estimate_xpts_sfm(eqs{k}, sims{k});    
    if isnan(xps{k+1}(1)), xps{k+1} = xps{k} + [-.005 .005 0 0]; end    
    eqs{k+1}  = designeq_ml(xps{k+1}, shot, time_ms, eqs{2});
    
  elseif sfp && ~constrain_sp
    xps{k+1} = estimate_xpts_sfp(eqs{k}, sims{k});
    eqs{k+1}  = designeq_ml(xps{k+1}, shot, time_ms, eqs{2});
    
  elseif sfp && constrain_sp
    sps = estimate_strike_pts(eqs{k}, sims{k}); 
    opts.constrain_sp = 1;
    opts.sp = sps;
    eqs{k+1} = designeq_ml(cake_xp, shot, time_ms, eqs{2}, opts);
    snf = analyzeSnowflake(eqs{k+1});
    xps{k+1} = [snf.rx snf.zx];                
  end
  
  sims{k+1} = heatsim_ml(eqs{k+1},shot,time_ms); 
end

eqs = {eqs{1}, eqs{2}, eqs{end}};

save('xps','xps')
save('eqs','eqs')
save('sims','sims')
































