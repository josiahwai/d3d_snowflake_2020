clear; clc; close all;

plotit = 1;
shot = 155355;
time_ms = 4400;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath('/u/jwai/d3d_snowflake_2020/current')); 
addpath(genpath('/u/jwai/d3d_snowflake_2020/test'));
load('/u/jwai/d3d_snowflake_2020/current/inputs/tok_data/d3d_obj_mks_struct_6565.mat')
efitdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/efit01/155355';
% hfdir = '/u/jwai/d3d_snowflake_2020/current/outputs/hfsims/efit_unconstrained/165288/';
% hffn = 'hfsim_165288_2900.mat';
% load([hfdir hffn]);
load('hfsim_155355_4400.mat')
eq = read_eq(shot,time_ms,efitdir);
eq = eq.gdata;

% ---------------------
% ANALYZE THE SNOWFLAKE
% ---------------------
limdata = tok_data_struct.limdata;

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

% xp = [[.05 -.05  .05 -.05] [ 1.0160  1.0364  1.1341  1.2653]-R0];
% vp = [[.05 -.05 -.05  .05] [-1.1281 -1.2462 -1.3439 -1.3643]-Z0];

xp = linspace(-.1,.1,5);
vp = linspace(-.1,.1,5);
c = combvec(xp,vp);
xp = c(1,:); vp = c(2,:);

np = length(xp);

if plotit
%   plot_eq(eq)
  plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
  hold on
  axis equal
  axis([1.0 1.5 -1.4 -0.9])
   scatter(R0+xp, Z0+vp,'g')
end

% curl-free constraint
Aeq1 = ...
  [0 -1 0 2*R0 0 2*R0 0 0 0 0; 
  0 0 0 0 0 2*R0 6*R0^2 0 2*R0^2 0;
  0 0 0 0 -2*R0 0 0 2*R0^2 0 6*R0^2];
  
beq1 = [0 0 0]';

% x-point constraint
xx = [rxPL rxPL rxSL rxSL] - R0;
vx = [zxPL zxPL zxSL zxSL] - Z0;

Aeq2 = gradpolyterms(xx,vx);
beq2 = [0 0 0 0]';


% probe the plasma at np places
[dpdx, dpdv] = gradpolyterms(xp,vp);
p = polyterms(xp,vp);

[psip, psip_r, psip_z] = bicubicHermite(rg,zg,psizr,R0+xp, Z0+vp);

Ap = [p; dpdx; dpdv];
bp = [psip psip_r psip_z]';


% find the coeff
Aeq = [Aeq1; Aeq2]; 
beq = [beq1; beq2];
c = leastsquare(Ap,bp,Aeq1,beq1);


% fitted flux (psizr2)
xg = rg - R0;
vg = zg - Z0;
psizr2 = zeros(length(vg),length(xg));

for j = 1:length(vg)
  for i = 1:length(xg)
    p = polyterms(xg(i),vg(j));
    psizr2(j,i) = p*c;
  end
end

% analyze the fitted flux
% .......................

% get a good starting pt for snowfinder
[r0,z0] = isoflux_xpFinder(psizr2,1.15,-1.25,rg,zg); 
 
[rxPL2, zxPL2, rxSL2, zxSL2] = snowFinder(psizr2,r0,z0,0.1,rg,zg); 

% zoom in on snowflake x-pts
[rxPL2, zxPL2, psixPL2] = isoflux_xpFinder(psizr2, rxPL2, zxPL2, rg, zg);
[rxSL2, zxSL2, psixSL2] = isoflux_xpFinder(psizr2, rxSL2, zxSL2, rg, zg);

contour(rg, zg, psizr2, [psixPL2 psixPL2], 'b')
contour(rg, zg, psizr2, [psixSL2 psixSL2], 'b')



% ---------------------------------------
% ESTIMATE NEW FLUX DISTRIBUTION FROM IR
% ---------------------------------------

% where is the peak based on IR
% .............................
qperp_dir  = [root 'inputs/qperp/' num2str(shot) '/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t

[~,k] = min(abs(t-time_ms));
qir = qperp(k,:)'/100;
s = s/100;

% Remove the gap from s (distance along limiter)
iGap = find(s < 1.70,1,'last');
dgap = s(iGap+1) - s(iGap);
s(iGap(end)+1:end) = s(iGap(end)+1:end) - dgap;

pkthresh = 1.5*median(qir);

iI = find(s<1.15);
iO = find(s>1.45);
iX = setdiff(1:length(s), [iI iO]);

sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

[~,s_qirmaxI] = qpeak_info(...
    qir(iI), s(iI), pkthresh, sLimTot, limdata, rg, zg, psizr2);

[~,s_qirmaxX] = qpeak_info(...
    qir(iX), s(iX), pkthresh, sLimTot, limdata, rg, zg, psizr2);

[~,s_qirmaxO] = qpeak_info(...
    qir(iO), s(iO), pkthresh, sLimTot, limdata, rg, zg, psizr2);

  
% where is the peak based on hf sim?
% .................................

[~,k] = max(hfsim.qdiv_perpI);
s_qmaxI = sLimTot - calcLimDistance(hfsim.rdivI(k), hfsim.zdivI(k), limdata);

[~,k] = max(hfsim.qdiv_perpX);
s_qmaxX = sLimTot - calcLimDistance(hfsim.rdivX(k), hfsim.zdivX(k), limdata);

[~,k] = max(hfsim.qdiv_perpO);
s_qmaxO = sLimTot - calcLimDistance(hfsim.rdivO(k), hfsim.zdivO(k), limdata);

% how much should the strike points be moved
ds = [s_qirmaxI s_qirmaxX s_qirmaxX s_qirmaxO] - [s_qmaxI s_qmaxX s_qmaxX s_qmaxO];
ds(4) = ds(4);
ds(1) = ds(1);

Wsp = diag([2 2 20 20]);

% find the current strike points 
% ..............................

% isoflux_spFinder is not robust and can only find it you start super close
[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 500, true);
psilim = bicubicHermite(rg,zg,psizr2,rlim,zlim);
[~,k] = findpeaks(-abs(psilim - psixSL2));

k = k(rlim(k) < 1.36 & rlim(k) > 0.8 & zlim(k) < -1 & zlim(k) > -1.4);

sprzI  = isoflux_spFinder(psizr2,psixPL2,rg,zg,[zlim rlim]',k(4)-5:k(4)+5);
sprzXI = isoflux_spFinder(psizr2,psixPL2,rg,zg,[zlim rlim]',k(3)-5:k(3)+5);
sprzXO = isoflux_spFinder(psizr2,psixSL2,rg,zg,[zlim rlim]',k(2)-5:k(2)+5);
sprzO  = isoflux_spFinder(psizr2,psixSL2,rg,zg,[zlim rlim]',k(1)-5:k(1)+5);

s_spI  = sLimTot - calcLimDistance(sprzI(1), sprzI(2), limdata);
s_spXI = sLimTot - calcLimDistance(sprzXI(1), sprzXI(2), limdata);
s_spXO = sLimTot - calcLimDistance(sprzXO(1), sprzXO(2), limdata);
s_spO  = sLimTot - calcLimDistance(sprzO(1), sprzO(2), limdata);

s_sp = [s_spI s_spXI s_spXO s_spO];


% desired strike points
s_spd = s_sp + ds;
r_spd = []; z_spd = [];
for k = 1:4
  [r_spd(k), z_spd(k)] = calcLimDistanceInv(sLimTot - s_spd(k), limdata);
end

scatter(r_spd,z_spd,'k','filled')

% strike point polynomials
p_spd = polyterms(r_spd - R0, z_spd - Z0);


for iter = 1:10

  % fitted x-points in linearized coordinate system
  xx = [rxPL2 rxPL2 rxSL2 rxSL2] - R0;
  vx = [zxPL2 zxPL2 zxSL2 zxSL2] - Z0;
  p_xp = polyterms(xx,vx);
  
  
  
  % Constrain
  
  Wc  = diag([1 R0 R0 R0^2 R0^2 R0^2 R0^3 R0^3 R0^3 R0^3]);
%   Wc = eye(10);
  
  Air = blkdiag(Wsp,Wc)*[p_spd - p_xp; eye(10)];
  bir = [zeros(4,1); Wc*c];
  
  c_new = leastsquare(Air,bir,Aeq1,beq1);

  % -------------------------------
  % ANALYZE THE IR-CONSTRAINED FLUX
  % -------------------------------
  
  for j = 1:length(vg)
    for i = 1:length(xg)
      p = polyterms(xg(i),vg(j));
      psizr2(j,i) = p*c_new;
    end
  end
  
  % analyze the fitted flux
  % .......................
  
  [r0,z0] = isoflux_xpFinder(psizr2,1.15,-1.25,rg,zg);
  [rxPL2, zxPL2, rxSL2, zxSL2] = snowFinder(psizr2,r0,z0,0.1,rg,zg);
  
  % zoom in on snowflake x-pts
  [rxPL2, zxPL2, psixPL2] = isoflux_xpFinder(psizr2, rxPL2, zxPL2, rg, zg);
  [rxSL2, zxSL2, psixSL2] = isoflux_xpFinder(psizr2, rxSL2, zxSL2, rg, zg);
  
end
% scatter(rxPL2, zxPL2,'r','filled')
% scatter(rxSL2, zxSL2,'r','filled')
contour(rg, zg, psizr2, [psixPL2 psixPL2], 'r')
contour(rg, zg, psizr2, [psixSL2 psixSL2], 'r')


% ---------------------------
% stuff
% ---------------------------

load('/u/jwai/d3d_snowflake_2020/current/optimize/convergence/interior-point/history.mat')

x0 = [1.1056 1.1603 -1.2072 -1.2597];
[r(1) r(2) z(1) z(2)] = unpack(x0);
plot(r, z, 'xb', 'Markersize', 10, 'LineWidth', 4)


[r(1) r(2) z(1) z(2)] = unpack(history.x(end,:));
plot(r, z, 'xk', 'Markersize', 10, 'LineWidth', 4)

plot([rxPL2 rxSL2], [zxPL2 zxSL2], 'xr', 'Markersize', 10, 'LineWidth', 4)

set(gcf,'position',[1300 383 560 420])

xp = [rxPL rxSL zxPL zxSL];















