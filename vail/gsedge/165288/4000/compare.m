
load eq_before.mat

jpar_before = eq_before.jpar;

psi1 = eq_before.psibar;
load eq_after.mat

jpar_after = eq_after.jpar;

psi2 = eq_after.psibar;

plot(psi1, jpar_before, '-r', 'LineWidth', 2)


hold on

plot(psi2, jpar_after, '-b', 'LineWidth', 2)

axis tight

xlabel('psi_N')
ylabel('J')

shot = 165288;
time = 4000;

% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

limdata = limdata(:,[1:79 81:end]);

% Define strike point segments
 
Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

zShelf = limdata(1,limIdxL(3));

% Compute total length of the limiter [m]

sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% Compute distance to limiter landmarks [m]

s45Deg1 = sLimTot - calcLimDistance(limdata(2,79), limdata(1,79), limdata);
s45Deg2 = sLimTot - calcLimDistance(limdata(2,78), limdata(1,78), limdata);

%...............................................................................
% Configure the plots

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
title([int2str(shot) ': ' int2str(time) ' ms'])


% Plot the initial equilibrium

efit_dirname = ['/u/pvail/d3d_snowflake_2019/gsedge/165288/4000' ];
    
eq = read_eq(shot, time/1000, efit_dirname);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
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




% Plot the final equilibrium

load eq_after.mat;

rg = eq_after.rg; 
zg = eq_after.zg;

psizr  = eq_after.psizr;
psibry = eq_after.psibry;

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


contour(rg, zg, psizr, [psixPL psixPL], 'r', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'r', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xr', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xr', 'Markersize', 12, 'LineWidth', 3)

figure(2)
%plot(jpar



