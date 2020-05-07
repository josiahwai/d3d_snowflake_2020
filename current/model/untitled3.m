% ========
% Settings
% ========
close all
shot = 155354;
time_ms = 4823;
xp1 = [1.1506    1.1856   -1.1631   -1.3263];
root = '/u/jwai/d3d_snowflake_2020/current/';

% =============
% Simulate eq0
% ============

% load
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
eq0 = read_eq(shot, time_ms/1000, cake_dir);
eq0 = eq0.gdata;
cake_snow = analyzeSnowflake(eq0);
xp0 = [cake_snow.rx cake_snow.zx];

% simulate
sim0 = heatsim_ml(eq0,shot,time_ms,1);
J0 = measure_cost4(sim0);

% ==========
% Find dxpP
% ==========
struct_to_ws(eq0); clear xlim ylim
struct_to_ws(sim0);

[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(1),z_qmax(1));
df1dr = -dpsidr / dpsidz; 

[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(2),z_qmax(2));
df2dr = -dpsidr / dpsidz; 


dz = [z_qirmax(1:2) - z_qmax(1:2)]';  % why doesn't this method work?
dr = [r_qirmax(1:2) - r_qmax(1:2)]';
dfdr = [df1dr; df2dr];

relax = .4;
b =  [dz - dfdr.*dr; 0; 0];
A = [-df1dr 1; -df2dr 1; relax*eye(2)];

erz = pinv(A)*b

rx1 = xp0(1) + erz(1);
zx1 = xp0(3) + erz(2);

scatter(rx1,zx1,'filled')











