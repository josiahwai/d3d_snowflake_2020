function plot_eq(eq,fignum)

if nargin < 2, fignum = 18; end

figure(fignum)
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
[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);

[~,~,~,~, psixPL, psixSL] = my_snowfinder(rg, zg, psizr, eq.psibry);

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
contour(rg,zg,psizr,[psixPL psixPL], 'b', 'linewidth', 1);
contour(rg,zg,psizr,[psixSL psixSL], 'b', 'linewidth', 1);

axis equal
% axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
    
    
    
    
    
    
    
    