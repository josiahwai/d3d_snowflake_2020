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

u = [dpsidz -dpsidr]' / norm([dpsidz -dpsidr]);


dz = [z_qirmax(1:2) - z_qmax(1:2)]';  % why doesn't this method work?
dr = [r_qirmax(1:2) - r_qmax(1:2)]';
dfdr = [df1dr; df2dr];

erz = inv([df1dr 1; df2dr 1]) * [dz - dfdr.*dr]

rprobe = rx(1) + [-.1:.001:.1];
zprobe = (zx(1) + dz(1)) * ones(1,length(rprobe));
psi = bicubicHermite(rg,zg,psizr,rprobe,zprobe);
[~,k] = min(abs(psi - psix(1)));
scatter(rprobe(k),zprobe(k),'filled')


v = [dr(2) dz(2)]';
r = v - u*u'*v;

xp1 = [rprobe(k) zprobe(k)] + norm(r) * [1 0]

er = xp1(1) - xp0(1)
ez = xp1(2) - xp0(3)

% dr = .01;
% r1 = r_qmax(1) + dr;
% z1 = z_qmax(1) + df1dr*dr; 
% contour(rg,zg,psizr,[psi_qmax(1) psi_qmax(1)], 'k')
% 
% dr = -.1:.005:.1;
% r2 = r_qmax(2) + dr;
% z2 = z_qmax(2) + df2dr*dr; 
% scatter(r2,z2)
% contour(rg,zg,psizr,[psi_qmax(2) psi_qmax(2)], 'k')









% Iterate to find best xp
% =======================
% drxP = .02;
% dzxP = .04;
% drxS = -.02;
% dzxS = 0;
% 
% close all
% xp1 = xp0 + [drxP drxS dzxP dzxS];
% eq1 = designeq_ml(xp1,shot,time_ms);
% sim1 = heatsim_ml(eq1,shot,time_ms,1);
% J1 = measure_cost4(sim1);
% 
% 
% % Iterate to find best xp
% % =======================
% drxP = .002;
% dzxP = .012;
% drxS = -.02;
% dzxS = 0;
% 
% close all
% xp1 = xp1 + [drxP drxS dzxP dzxS];
% eq1 = designeq_ml(xp1,shot,time_ms);
% sim1 = heatsim_ml(eq1,shot,time_ms,1);
% J1 = measure_cost4(sim1);

