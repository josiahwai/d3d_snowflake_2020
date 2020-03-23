clear; clc; close all;

plotit = 1;

addpath(genpath('/u/jwai/d3d_snowflake_2020/current')); 
addpath(genpath('/u/jwai/d3d_snowflake_2020/test'));
load('/u/jwai/d3d_snowflake_2020/current/inputs/tok_data/d3d_obj_mks_struct_6565.mat')
efitdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/efit01/165288';

eq = read_eq(165288,4000,efitdir);
eq = eq.gdata;


% ---------------------
% ANALYZE THE SNOWFLAKE
% ---------------------

% Locate the snowflake in the lower divertor
% ..........................................
psibry = eq.psibry;
[psizr, rg, zg] = regrid(eq.rg, eq.zg, eq.psizr, 257, 257);

% get a good starting pt for snowfinder
[r0,z0] = isoflux_xpFinder(psizr,1.15,-1.25,rg,zg); 

[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, r0, z0, 0.1, rg, zg); 

% zoom in on snowflake x-pts
[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    swap(psixPL, psixSL);
    swap(rxPL, rxSL);
    swap(zxPL, zxSL);   
end

% ---------------------------------
% 3RD ORDER POLYNOMIAL FIT TO psizr
% ---------------------------------

% polynomial expansion of flux
R0 = (rxPL + rxSL)/2;
Z0 = (zxPL + zxSL)/2;

xp = [.05 -.05  .05 0];
vp = [.05 -.05 -.05 0];
np = length(xp);

if plotit
  plot_eq(eq)
  axis([1.0 1.5 -1.4 -0.9])
  scatter(R0+xp, Z0+vp,'k','filled')
end

% const, l1, l2, q1, q2, q3, c1, c2, c3, c4
Aeq = ...
  [0 -1 0 2*R0 0 2*R0 0 0 0 0; 
  0 0 0 0 0 2*R0 6*R0^2 0 2*R0^2 0;
  0 0 0 0 -2*R0 0 0 2*R0^2 0 6*R0^2];
  
beq = [0 0 0]';

% probe the plasma at 3 places
[dpdx, dpdv] = gradpolyterms(xp,vp);
p = polyterms(xp,vp);

[psip, psip_r, psip_z] = bicubicHermite(rg,zg,psizr,R0+xp, Z0+vp);

Ap = [p; dpdx; dpdv];
bp = [psip psip_r psip_z]';


% minimize over c, f(c) = norm(Ap*c-bp) subject to Aeq*c = beq
% c := 10 taylor expansion coeffiecients  = const, l1, l2 ... c3, c4

dum = pinv([Ap'*Ap Aeq'; Aeq zeros(3,3)]) * [Ap'*bp; beq];
c = dum(1:10);  
lam = dum(11:end); % lagrange multipliers


% fitted flux (psizr2)
xg = rg(rg >  1   & rg <  1.5) - R0;
vg = zg(zg > -1.4 & zg < -0.9) - Z0;

psizr2 = zeros(length(vg),length(xg));

for j = 1:length(vg)
  for i = 1:length(xg)
    p = polyterms(xg(i),vg(j));
    psizr2(j,i) = p*c;
  end
end

% analyze the fitted flux
% .......................
rg = xg + R0;
zg = vg + Z0;

% get a good starting pt for snowfinder
[r0,z0] = isoflux_xpFinder(psizr,1.15,-1.25,rg,zg); 

[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr2,r0,z0,0.1,rg,zg); 

% zoom in on snowflake x-pts
[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr2, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr2, rxSL, zxSL, rg, zg);

contour(xg+R0, vg+Z0, psizr2, [psixPL psixPL], 'k','linewidth',2)
contour(xg+R0, vg+Z0, psizr2, [psixSL psixSL], 'k')

% contour(xg+R0, vg+Z0, psizr2, 200)




































