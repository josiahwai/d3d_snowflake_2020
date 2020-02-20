
time = 4000;

shot = 165288;

%...............................................................................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

limdata = limdata(:,[1:79 81:end]);

%...............................................................................
% Configure the plots

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
% axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
title([int2str(shot) ': ' int2str(time) ' ms'])

%...............................................................................
% Load and analyze the equilibrium

eqdir = '/u/pvail/d3d_snowflake_2019/heatfluxconstrain/165288/4000/';
eqnam = 'eq_165288_4000_EQHF.mat';

load([eqdir eqnam])

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

contour(rg, zg, psizr, [psixPL psixPL], 'k', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'k', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xk', 'Markersize', 12, 'LineWidth', 3)

%...............................................................................
% Load and analyze the equilibrium

load('eq_62.mat')

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

contour(rg, zg, psizr, [psixPL psixPL], 'm', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'm', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xm', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xm', 'Markersize', 12, 'LineWidth', 3)

%...............................................................................
% Load and analyze the equilibrium

load('eq_63.mat')

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

contour(rg, zg, psizr, [psixPL psixPL], 'g', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'g', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xg', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xg', 'Markersize', 12, 'LineWidth', 3)

%...............................................................................
% Load and analyze the equilibrium

load('eq_61.mat')

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

contour(rg, zg, psizr, [psixPL psixPL], 'c', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'c', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xc', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xc', 'Markersize', 12, 'LineWidth', 3)

