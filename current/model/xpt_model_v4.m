close all
root = '/u/jwai/d3d_snowflake_2020/current/';
shot = 155354;
time_ms = 4823;


% =============
% Simulate eq0
% ============

% load and simulate
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
eq0 = read_eq(shot, time_ms/1000, cake_dir);

eq0 = eq0.gdata;

cake_snow = analyzeSnowflake(eq0);
xp0 = [cake_snow.rx cake_snow.zx];
% eq0 = designeq_ml(xp0,shot,time_ms);

close all
sim0 = heatsim_ml(eq0,shot,time_ms,1);
close(11)


% =========================
% Find new x-pts & simulate
% =========================
close all
xp1 = estimate_xpts(eq0,sim0,1)
eq1 = designeq_ml(xp1,shot,time_ms);
sim1 = heatsim_ml(eq1,shot,time_ms,1);

close all
xp2 = estimate_xpts(eq1,sim1,1)
eq2 = designeq_ml(xp2,shot,time_ms);
sim2 = heatsim_ml(eq2,shot,time_ms,1);


close all
xp3 = estimate_xpts(eq2,sim2,1)
eq3 = designeq_ml(xp3,shot,time_ms);
sim3 = heatsim_ml(eq3,shot,time_ms,1);


close all
xp4 = estimate_xpts(eq3,sim3,1)
eq4 = designeq_ml(xp4,shot,time_ms);
sim4 = heatsim_ml(eq4,shot,time_ms,1);






















