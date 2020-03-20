% snowfinder 2.0
ccc
load('eq165288_2900_v7.mat')

% ---------
% PLOT INFO
% ---------

% Load tokamak definition
root = '/u/jwai/d3d_snowflake_2020/current/';
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);

limdata = tok_data_struct.limdata;
[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 500, true);

figure(11)
plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([0.9 1.6 -1.5 -0.8])
xlabel('R [m]')
ylabel('Z [m]')

%------------------
% SNOWFLAKE TESTING
% -----------------
psibry = eq.psibry;
[psizr,rg, zg] = regrid(eq.rg, eq.zg, eq.psizr, 257, 257);

[rExp,zExp] = isoflux_xpFinder(psizr,1.15,-1.25,rg,zg);

% rExp   =  1.1300;
% zExp   = -1.2200;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);
psixPL = bicubicHermite(rg,zg,psizr,rxPL,zxPL);
psixSL = bicubicHermite(rg,zg,psizr,rxSL,zxSL);




[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

% scatter(rExp,zExp,'k','filled')
% scatter(rxPL,zxPL,'r','filled')
% scatter(rxSL,zxSL,'b','filled')

scatter(rxPL,zxPL,'r','filled')
scatter(rxSL,zxSL,'b','filled')
contour(rg, zg, psizr, [psixPL psixPL], 'r', 'LineWidth', 2)
contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)



































