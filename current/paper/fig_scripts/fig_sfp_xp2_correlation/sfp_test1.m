clear all; clc; close all; warning('off','all');

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% =========================
% LOAD AND SIMULATE EFIT EQS
% =========================

efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

% simulate the efit01
efit_eq1 = read_eq(shot, time_ms/1000, efit_dir);
efit_snow1 = analyzeSnowflake(efit_eq1);
efit_xp1 = [efit_snow1.rx efit_snow1.zx];

% simulate a hybrid eq: a toksys eq with x-pts defined by efit eq
clear iter init spec config gsdesign 
efit_xp2 = efit_xp1;
efit_eq2 = designeq_ml(efit_xp2,shot,time_ms);
eq0 = efit_eq2;

% CORRELATION TESTS
N = 16;
th = linspace(0,2*pi,N+1);
th(end) = [];

snow0 = analyzeSnowflake(eq0);
xp0 = [snow0.rx snow0.zx];
sp0 = [snow0.rSPP snow0.zSPP];

sps = []; 
xps = [];
eqs = {};

for k = 1:2*N
  
  if k < N
    dxp(k,:) = (0.5 - k / N) * .04 * [0 1 0 0];
  else
    dxp(k,:) = (1.5 - k / N) * .04 * [0 0 0 1];
  end
  
  xps(k,:) = xp0 + dxp(k,:);
  
  eqs{k} = designeq_ml( xps(k,:), shot, time_ms);
  
  snow = analyzeSnowflake(eqs{k});
  
  sps(k,:) = [snow.rSPP snow.zSPP];
  
  dsp(k,:) = sps(k,:) - sp0;
end


% save it
savedir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/sfp_test/';

save([savedir + 'dxp'], 'dxp')
save([savedir + 'dsp'], 'dsp')




















