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

%==========================================================================
%==========================================================================

load('/u/jwai/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKE/eq_165288_4000_CAKE.mat')

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
contour(rg, zg, psizr, [psixPL psixPL], '-b', ...
    'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], '-b', ...
    'LineWidth', 1)

plot(rxPL, zxPL, 'xb', 'Markersize', 12, ...
    'LineWidth', 3)
plot(rxSL, zxSL, 'xb', 'Markersize', 12, ...
    'LineWidth', 3)


%==========================================================================
%==========================================================================




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

labels = {'Unconstrained kEFIT', 'Constrained kEFIT'};
co = {[0 0 1], [0.91 0.41 0.17]};
location = 'northeast';
mylegend(labels,co,location)













