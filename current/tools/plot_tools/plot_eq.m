
if isfield(eq,'gdata')
    eq = eq.gdata;
end

root = '/u/jwai/d3d_snowflake_2020/current/';
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);

limdata = tok_data_struct.limdata;
psizr = eq.psizr;
rg = eq.rg';
zg = eq.zg;
rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
contour(rg,zg,psizr,[psixPL psixPL], 'r', 'linewidth', 2);
contour(rg,zg,psizr,[psixSL psixSL], '--b', 'linewidth', 1);

axis equal
% axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
    
    
    
    
    
    
    
    