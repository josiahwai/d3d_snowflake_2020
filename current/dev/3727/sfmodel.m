clear all; clc; close all; warning('off','all');

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% 155328: .009, .004, 0.2, 0.05
% 155355: .006, .002, 0.07, 0.02

lambdaq_i = .0063;
lambdaq_o = .0023;
chi_i = .1;
chi_o = .03;

% =========================
% LOAD AND SIMULATE EFIT EQS
% =========================

% load and simulate
efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

% simulate the efit01
efit_eq1 = read_eq(shot, time_ms/1000, efit_dir);
efit_sim1 = 0; % heatsim_fit(efit_eq1.gdata, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);
efit_snow1 = analyzeSnowflake(efit_eq1);
efit_xp1 = [efit_snow1.rx efit_snow1.zx];
xp0 = efit_xp1;

% simulate a hybrid eq: a toksys eq with x-pts defined by efit eq
clear iter init spec config gsdesign 
efit_xp2 = efit_xp1;
efit_eq2 = designeq_ml_beta(efit_xp2,shot,time_ms);
efit_sim2 = heatsim_fit(efit_eq2, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);
efit_eq2 = rmfield(efit_eq2, {'r', 'b', 'p'});

% ============================
% Fit the Eich Profile to data
% ============================
load('d3d_obj_mks_struct_6565.mat')

% Load heat flux data q(s,t), s=distance along limiter, and t=time
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
[~,k] = min(abs(t-time_ms));
qperp = qperp(k,:)';

% obtain parameters from the eich fit to heat flux (strike points etc.)
ef = eich_fitter(s', qperp, efit_eq1, tok_data_struct);


% =========================
% Find new x-pts & simulate
% =========================
N = 11;
if constrain_sp, N = 2; end

eqs = {efit_eq1, efit_eq2};
sims = {efit_sim1, efit_sim2};
xps  = {efit_xp1, efit_xp2};

sfp = 0;
if isnan(ef.ssp(2)), sfp = 1; end
sfm = ~sfp;


for k = 2:N
  fprintf(['\nIteration ' num2str(k) ' of ' num2str(N)])
  
  % constrain x-pts in snowflake minus
  if sfm
    xps{k+1}  = estimate_xpts_sfm(eqs{k}, ef, sims{k});
    eqs{k+1}  = designeq_ml_beta(xps{k+1}, shot, time_ms, eqs{k});
    
  % constrain x-pts in snowflake plus
  elseif sfp && ~constrain_sp
    xps{k+1} = estimate_xpts_sfp5(eqs{k}, ef);
    eqs{k+1}  = designeq_ml_beta(xps{k+1}, shot, time_ms, eqs{k});
    
  % constrain strike-pts in snowflake minus (1 iteration only)
  elseif sfp && constrain_sp    
    opts.constrain_sp = 1;
    opts.sp = [ef.rsp([1 end]) ef.zsp([1 end])];
    eqs{k+1} = designeq_ml_beta(xp0, shot, time_ms, eqs{k}, opts);
    
    snf = analyzeSnowflake(eqs{k+1});
    xps{k+1} = [snf.rx snf.zx];    
    sims{k+1} = heatsim_fit(eqs{k+1}, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);
    break
  end
  
  sims{k+1} = heatsim_fit(eqs{k+1}, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);
 
  % clear memory
  eqs{k+1} = rmfield(eqs{k+1}, {'r', 'b', 'p'});
end

eqs(3:end-1) = [];
xps(3:end-1) = [];
sims(3:end-1) = [];

save('xps','xps')
save('sims','sims')
save('eqs', 'eqs')

delete('*.out')
delete('*.err')






























