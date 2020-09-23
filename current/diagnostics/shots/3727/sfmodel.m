clear all; clc; close all; warning('off','all');

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% 155328: .009, .004, 0.2, 0.05
% 155355: .006, .002, 0.07, 0.02

lambdaq_i = .0063;
lambdaq_o = .0018;
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
% efit_sim1 = heatsim_fit(efit_eq1.gdata, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);
efit_snow1 = analyzeSnowflake(efit_eq1);
efit_xp1 = [efit_snow1.rx efit_snow1.zx];
xp0 = efit_xp1;

% simulate a hybrid eq: a toksys eq with x-pts defined by efit eq
clear iter init spec config gsdesign 
efit_xp2 = efit_xp1;
efit_eq2 = designeq_bp(efit_xp2,shot,time_ms);
% efit_sim2 = heatsim_fit(efit_eq2, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o);


load('xps')
xpf = xps{end};
eqf = designeq_bp(xpf,shot,time_ms);

% eq0 = efit_eq1;
% eq1 = rmfield(efit_eq2, {'r', 'b', 'p'});
eqf = rmfield(eqf, {'r', 'b', 'p'});

% save('eq0','eq0')
% save('eq1','eq1')
% save('eqf','eqf')

compare_bp_fl_eqs



























