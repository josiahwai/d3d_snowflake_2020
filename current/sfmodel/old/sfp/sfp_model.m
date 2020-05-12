clear all; clc; close all
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root))
shot = 155354;
time_ms = 2730;

% 2000 - reverse
% 3528 - transition
% 3627 - transition
% 4125 - close


% =============
% Simulate eq0
% ============

% load and simulate
try 
  load('155354_2730_eq0')
  load('155354_2730_sim0.mat')
catch  
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  eq0 = read_eq(shot, time_ms/1000, cake_dir);
  eq0 = eq0.gdata;    
  snow0 = analyzeSnowflake(eq0);
  xp0 = [snow0.rx snow0.zx];
  eq0 = designeq_ml(xp0,shot,time_ms);
  sim0 = heatsim_ml(eq0,shot,time_ms);
end

snow0 = analyzeSnowflake(eq0);
xp0 = [snow0.rx snow0.zx];



% =========================
% Find new x-pts & simulate
% =========================

N = 10; 

eqs{1}  = eq0;
sims{1} = sim0;
xps{1}  = xp0;

for k = 1:N
  k
  xps{k+1} = estimate_xpts_sfp(eqs{k},sims{k},1);
  eqs{k+1} = designeq_ml(xps{k+1},shot,time_ms,eqs{k});
  sims{k+1} = heatsim_ml(eqs{k+1},shot,time_ms); 
end

% =========
% Testing
% =========

ccc
load('155354_2730_eq0')
shot = 155354;
time_ms = 2730;
snow0 = analyzeSnowflake(eq0);
xp0 = [snow0.rx snow0.zx];

eq1 = designeq_ml(xp0, shot, time_ms, eq0);

load('155354_2730_eq0')
load('155354_2730_sim0.mat')
struct_to_ws(sim0);
snow0 = analyzeSnowflake(eq0);
rsppi = snow0.rSPP(1);
zsppi = snow0.zSPP(1) - s_qirmax(1) + s_qmax(1);
rsppo = snow0.rSPP(2) + s_qirmax(3) - s_qmax(3);
zsppo = snow0.zSPP(2);

close all
plot_eq(eq1)
axis([1 1.5 -1.4 -.9])
scatter([rsppi rsppo], [zsppi zsppo], 'filled')

sim1 = heatsim_ml(eq1,shot,time_ms);
plotsim(sim1)

struct_to_ws(eq1);
contour(rg,zg,psizr,[psibry psibry], 'k', 'linewidth', 2)

snow1 = analyzeSnowflake(eq1);
xp1 = [snow1.rx snow1.zx];
xp1 - xp0

load('xps')
xps{4} - xp0

xps{4} - xp1

% xp1 = xp0;
% dxp = estimate_xpts_sfp(eq0,sim0,plotit)
% xp1([1 3]) = xp0([1 3]) + dxp';
% eq1 = designeq_ml(xp1, shot, time_ms, eq0);
% sim1 = heatsim_ml(eq1,shot,time_ms);
% 
% xp2 = xp1;
% dxp = estimate_xpts_sfp(eq1,sim1)
% xp2([1 3]) = xp1([1 3]) + dxp';
% eq2 = designeq_ml(xp2, shot, time_ms, eq0);
% sim2 = heatsim_ml(eq2,shot,time_ms);
% 
% 
% % dxp = estimate_xpts_sfp(eq2,sim2)
% % xp3([1 3]) = xp2([1 3]) + dxp';
% 
% 
% xp2 = [1.2032    1.1269   -1.0981   -1.3786];
% xp3 = xp2;
% xp3([2 4]) = ginput
% 
% xp3 = [1.2032    1.1523   -1.0981   -1.3554];
% 
% eq3 = designeq_ml(xp3, shot, time_ms, eq0);
% sim3 = heatsim_ml(eq3,shot,time_ms);
% 
% dxp = estimate_xpts_sfp(eqs,sims)
% xps = xp3;
% xps([1 3]) = xp3([1 3]) + dxp';
% 
% xps([1 3]) = ginput
% xps([2 4]) = ginput
% 
% eqs = designeq_ml(xps,shot,time_ms,eq0);
% sims = heatsim_ml(eqs,shot,time_ms,1);
% plotsim(sims)
% 
% figure
% plot_eq(eqs)
% axis([1 1.5 -1.4 -.9])
% 
% scatter(xps(1:2),xps(3:4),'filled')
% 
% close all
% plotsim(sim0)
% plotsim(sim1)
% plotsim(sim2)
% 
% figure
% plot_eq(eq0)
% % plot_eq(eq1)
% axis([1 1.5 -1.4 -.9])
































