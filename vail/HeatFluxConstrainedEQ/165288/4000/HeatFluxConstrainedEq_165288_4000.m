
% Generate equilibrium for DIII-D with IRTV-derived constraint on SFD nulls

%........................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

% Compute total length of the limiter [m]

sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% Configure the plots

figure(11)
plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
title([int2str(165288) ': ' int2str(4000) ' ms'])
    
% Load the equilibrium

top_dir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/';

eq_dir = [top_dir 'EFIT01/'];

eq = read_eq(165288, 4.00, eq_dir);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

rmaxis = eq.gdata.rmaxis;
zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
psimag = eq.gdata.psimag;
psibry = eq.gdata.psibry;

[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);

% Locate the snowflake in the lower divertor

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    
    temp   = psixPL;
    psixPL = psixSL;
    psixSL = temp;
    
    temp1 = rxPL;
    temp2 = zxPL;
    rxPL  = rxSL;
    zxPL  = zxSL;
    rxSL  = temp1;
    zxSL  = temp2;
    
end

contour(rg, zg, psizr, [psixPL psixPL], 'b', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xb', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xb', 'Markersize', 12, 'LineWidth', 3)

% Find the strike points in the EFIT

Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);

r_SP1_EFIT01 = spRZLP(1,1); 
r_SP2_EFIT01 = spRZLP(2,1); 
r_SP3_EFIT01 = spRZLS(4,1);

z_SP1_EFIT01 = spRZLP(1,2); 
z_SP2_EFIT01 = spRZLP(2,2); 
z_SP3_EFIT01 = spRZLS(4,2);

plot(r_SP1_EFIT01, z_SP1_EFIT01, 'ob', 'LineWidth', 2)
plot(r_SP2_EFIT01, z_SP2_EFIT01, 'ob', 'LineWidth', 2)
plot(r_SP3_EFIT01, z_SP3_EFIT01, 'ob', 'LineWidth', 2)

fExp_SP1_EFIT01 = calcFluxExpansion(r_SP1_EFIT01, z_SP1_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);
fExp_SP2_EFIT01 = calcFluxExpansion(r_SP2_EFIT01, z_SP2_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);
fExp_SP3_EFIT01 = calcFluxExpansion(r_SP3_EFIT01, z_SP3_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);

s_SP1_EFIT01 = sLimTot - calcLimDistance(r_SP1_EFIT01, z_SP1_EFIT01, limdata); 
s_SP2_EFIT01 = sLimTot - calcLimDistance(r_SP2_EFIT01, z_SP2_EFIT01, limdata);
s_SP3_EFIT01 = sLimTot - calcLimDistance(r_SP3_EFIT01, z_SP3_EFIT01, limdata);

%...............................................
% Load the heat flux data for DIII-D shot 165288

datadir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/dataHF/';

datanam_q = 'qperp_165288_4000.txt';
datanam_s = 's_165288.txt';

qperp = transpose(importdata([datadir datanam_q]));
s = transpose(importdata([datadir datanam_s]));

% Remove the gap

idxGap = find(s < 170);

% Remove the gap

gap1 = s(idxGap(end));
gap2 = s(idxGap(end)+1);

dgap = gap2 - gap1;

s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;

% Index data for each SP

idx_SP1 = find(s < 115);
idx_SP3 = find(s > 145);

idx_SP2 = setdiff(1:length(s), [idx_SP1 idx_SP3]);

qperp_SP1 = qperp(idx_SP1);
qperp_SP2 = qperp(idx_SP2);
qperp_SP3 = qperp(idx_SP3);

s_SP1 = s(idx_SP1)';
s_SP2 = s(idx_SP2)';
s_SP3 = s(idx_SP3)';

figure(22)
plot(s, qperp, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
hold on
plot([50 210], [0 0], '--k')

axis([80 180 -2 50])
xlabel('s [cm]')
ylabel('Heat Flux W/cm^2')
title([int2str(165288) ': ' int2str(4000) ' ms'])

% Determine magnitude and location of heat flux peak at each SP

[qmax_SP1_temp, idx_SP1] = max(qperp_SP1);
[qmax_SP2_temp, idx_SP2] = max(qperp_SP2);
[qmax_SP3_temp, idx_SP3] = max(qperp_SP3);

sDiv_qmax_SP1_temp = s_SP1(idx_SP1);
sDiv_qmax_SP2_temp = s_SP2(idx_SP2);
sDiv_qmax_SP3_temp = s_SP3(idx_SP3);

% Determine FWHM at each SP

idx11 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'first');
idx21 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'last');

FWHM_SP1_temp = abs(s_SP1(idx11) - s_SP1(idx21));

idx12 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'first');
idx22 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'last');

FWHM_SP2_temp = abs(s_SP2(idx12) - s_SP2(idx22));

idx13 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'first');
idx23 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'last');

FWHM_SP3_temp = abs(s_SP3(idx13) - s_SP3(idx23));

sDiv_qmax_SP1_temp = 0.01*sDiv_qmax_SP1_temp;
sDiv_qmax_SP2_temp = 0.01*sDiv_qmax_SP2_temp;
sDiv_qmax_SP3_temp = 0.01*sDiv_qmax_SP3_temp;

FWHM_SP1_temp = 0.01*FWHM_SP1_temp;
FWHM_SP2_temp = 0.01*FWHM_SP2_temp;
FWHM_SP3_temp = 0.01*FWHM_SP3_temp;

qmax_SP1_temp = qmax_SP1_temp*(1e-6)*(100)^2;
qmax_SP2_temp = qmax_SP2_temp*(1e-6)*(100)^2;
qmax_SP3_temp = qmax_SP3_temp*(1e-6)*(100)^2;

%....................................................
% Fit Eich profile to IRTV data for each strike point

qEich_Inner = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (-(x-s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x-s0))./S) + qBG;

qEich_Outer = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;

qGauss = @(q0, s0, sigma, qBG, x) q0*exp(-(x-s0).^2/sigma) + qBG;

ftInner = fittype(qEich_Inner, 'problem', 'fExp', 'independent', 'x');
ftOuter = fittype(qEich_Outer, 'problem', 'fExp', 'independent', 'x');
ftGauss = fittype(qGauss,                         'independent', 'x');

optionsInner = fitoptions(ftInner);
optionsInner.Lower = [0 0 0 0 0];

optionsOuter = fitoptions(ftOuter);
optionsOuter.Lower = [0 0 0 0 0];

optionsGauss = fitoptions(ftGauss);
optionsGauss.Lower = [0 0 0 0];

%........................
% Eich profile fit to SP1

[maxq, idx_maxq] = max(qperp_SP1);
 
optionsInner.StartPoint = [2*maxq 1.4 0.2 s_SP1(idx_maxq) 0];

fitEich_SP1 = fit(s_SP1, qperp_SP1, ftInner, 'problem', fExp_SP1_EFIT01, ...
    optionsInner);
 
% Plot the optimized Eich profile

s_Eich = linspace(60,170,1000);
    
qperp_Eich_SP1 = eichProfileInnerWithBG(s_Eich, fitEich_SP1.q0,          ...
    fitEich_SP1.S, fitEich_SP1.lambdaQ, fExp_SP1_EFIT01, fitEich_SP1.s0, ...
    fitEich_SP1.qBG);
    
hold on
plot(s_Eich, qperp_Eich_SP1, '-r', 'LineWidth', 1)

plot([fitEich_SP1.s0 fitEich_SP1.s0], [-10 60], '--k', 'LineWidth', 1.5)

%.........................
% Gauss profile fit to SP2

[maxq, idx_maxq] = max(qperp_SP2);
 
optionsGauss.StartPoint = [2*maxq s_SP2(idx_maxq) 1 0]; 

fitGauss_SP2 = fit(s_SP2, qperp_SP2, ftGauss, optionsGauss);
 
% Plot the optimized Eich profile

s_Gauss = linspace(60,170,1000);
    
qperp_Gauss_SP2 = gaussProfileWithBG(s_Eich, fitGauss_SP2.q0, ...
    fitGauss_SP2.s0, fitGauss_SP2.sigma, fitGauss_SP2.qBG);

hold on
plot(s_Gauss, qperp_Gauss_SP2, '-r', 'LineWidth', 1)

plot([fitGauss_SP2.s0 fitGauss_SP2.s0], [-10 60], '--k', 'LineWidth', 1.5)

%........................
% Eich profile fit to SP3

[maxq, idx_maxq] = max(qperp_SP3);
 
optionsOuter.StartPoint = [2*maxq 1.4 0.2 s_SP3(idx_maxq) 0];

fitEich_SP3 = fit(s_SP3, qperp_SP3, ftOuter, 'problem', fExp_SP3_EFIT01, ...
    optionsOuter);
 
% Plot the optimized Eich profile

s_Eich = linspace(60,200,1000);
    
qperp_Eich_SP3 = eichProfileOuterWithBG(s_Eich, fitEich_SP3.q0,          ...
    fitEich_SP3.S, fitEich_SP3.lambdaQ, fExp_SP3_EFIT01, fitEich_SP3.s0, ...
    fitEich_SP3.qBG);
    
hold on
plot(s_Eich, qperp_Eich_SP3, '-r', 'LineWidth', 1)

plot([fitEich_SP3.s0 fitEich_SP3.s0], [-10 60], '--k', 'LineWidth', 1.5)

%..........................................
% Plot the strike point locations from IRTV

s_SP1_IRTV = fitEich_SP1 .s0;
s_SP2_IRTV = fitGauss_SP2.s0;
s_SP3_IRTV = fitEich_SP3 .s0;

sInv1 = sLimTot - 0.01*s_SP1_IRTV;
sInv2 = sLimTot - 0.01*s_SP2_IRTV;
sInv3 = sLimTot - 0.01*s_SP3_IRTV;

[r_SP1_IRTV, z_SP1_IRTV] = calcLimDistanceInv(sInv1, limdata);
[r_SP2_IRTV, z_SP2_IRTV] = calcLimDistanceInv(sInv2, limdata);
[r_SP3_IRTV, z_SP3_IRTV] = calcLimDistanceInv(sInv3, limdata);

figure(11)
plot(r_SP1_IRTV, z_SP1_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2)
plot(r_SP2_IRTV, z_SP2_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2)
plot(r_SP3_IRTV, z_SP3_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2)

%............................................
% Plot the strike point locations from EFIT01

figure(22)
plot(100*[s_SP1_EFIT01 s_SP1_EFIT01], [-10 60], '--b', 'LineWidth', 1.5)
plot(100*[s_SP2_EFIT01 s_SP2_EFIT01], [-10 60], '--b', 'LineWidth', 1.5)
plot(100*[s_SP3_EFIT01 s_SP3_EFIT01], [-10 60], '--b', 'LineWidth', 1.5)

%....................................................
% Predict the null locations using the heat flux data

rxP_PRED_TARGET =  1.1302 - 0.01;
zxP_PRED_TARGET = -1.1142 + 0.05;

rxS_PRED_TARGET =  1.1838 - 0.04;
zxS_PRED_TARGET = -1.2996 + 0.02;


rxP_PRED_TARGET =  1;
zxP_PRED_TARGET = -1;
rxS_PRED_TARGET =  1;
zxS_PRED_TARGET = -1;


figure(11)
plot(rxP_PRED_TARGET, zxP_PRED_TARGET, 'xg', 'Markersize', 12, 'LineWidth', 3)
hold on
plot(rxS_PRED_TARGET, zxS_PRED_TARGET, 'xg', 'Markersize', 12, 'LineWidth', 3)

HeatFluxMapping_RegressionTree_v6

% Scale the data to match the total power used in heat flux simulations

rxP_PRED = predict(rtree_rxP, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp])
rxS_PRED = predict(rtree_rxS, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp])
zxP_PRED = predict(rtree_zxP, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp])
zxS_PRED = predict(rtree_zxS, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp])

figure(11)
plot(rxP_PRED, zxP_PRED, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
    'LineWidth', 3)
hold on
plot(rxS_PRED, zxS_PRED, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
    'LineWidth', 3)

%..................................................................
% Design a Grad-Shafranov equilibrium with constraints on the nulls

addpath(genpath('/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/gatools'))

%gsdesign_d3d_165288_4000_EQHF
load('/u/jwai/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKEHF/eq_165288_4000_CAKEHF.mat')
rg = eq.rg; 
zg = eq.zg;

psizr  = eq.psizr;
psibry = eq.psibry;

[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);

% Locate the snowflake in the lower divertor

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    
    temp   = psixPL;
    psixPL = psixSL;
    psixSL = temp;
    
    temp1 = rxPL;
    temp2 = zxPL;
    rxPL  = rxSL;
    zxPL  = zxSL;
    rxSL  = temp1;
    zxSL  = temp2;
    
end

figure(11)
hold on
contour(rg, zg, psizr, [psixPL psixPL], 'color', [0.91 0.41 0.17], ...
    'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'color', [0.91 0.41 0.17], ...
    'LineWidth', 1)

plot(rxPL, zxPL, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
    'LineWidth', 3)
plot(rxSL, zxSL, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
    'LineWidth', 3)

% Find the strike points in the kinetic EFIT

Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);

r_SP1 = spRZLP(1,1); r_SP2 = spRZLP(2,1); r_SP3 = spRZLS(4,1);
z_SP1 = spRZLP(1,2); z_SP2 = spRZLP(2,2); z_SP3 = spRZLS(4,2);

% Compute distance to strike points [m]

s_SP1_EQHF = sLimTot - calcLimDistance(r_SP1, z_SP1, limdata); 
s_SP2_EQHF = sLimTot - calcLimDistance(r_SP2, z_SP2, limdata);
s_SP3_EQHF = sLimTot - calcLimDistance(r_SP3, z_SP3, limdata);

figure(22)
plot(100*[s_SP1_EQHF s_SP1_EQHF], [-10 60], '--', 'color', [0.91 0.41 0.17], ...
    'LineWidth', 1.5 )
plot(100*[s_SP2_EQHF s_SP2_EQHF], [-10 60], '--', 'color', [0.91 0.41 0.17], ...
    'LineWidth', 1.5 )
plot(100*[s_SP3_EQHF s_SP3_EQHF], [-10 60], '--', 'color', [0.91 0.41 0.17], ...
    'LineWidth', 1.5)
