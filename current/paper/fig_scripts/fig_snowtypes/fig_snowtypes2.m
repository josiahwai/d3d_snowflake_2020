close all
saveit = 0;


load('d3d_obj_mks_struct_6565.mat');
struct_to_ws(tok_data_struct);

dgray = [1 1 1]*0.4;
lgray = [1 1 1]*0.75;
limgray = [1 1 1] * 0.7;

% PLOT SETTINGS
ax(1) = subplot(221);
ax(2) = subplot(222);
ax(3) = subplot(223);
ax(4) = subplot(224);

ax(1).Position = [0.1 0.5 0.4 0.4];
ax(2).Position = [0.5 0.5 0.4 0.4];
ax(3).Position = [0.1 0.1 0.4 0.4];
ax(4).Position = [0.5 0.1 0.4 0.4];


set(ax, 'XTick', [], 'YTick', [], 'box', 'on')
set(ax, 'XLim', [0.8685    1.6246])
set(ax, 'YLim', [-1.4644  -0.8616])

set(gcf, 'position', [552 401 570 479])


% STANDARD
axes(ax(1))
hold on
load('eq_std')
struct_to_ws(eq.gdata);
[psizr, rg, zg] = regrid(rg, zg, psizr, 300, 300);

plot(limdata(2,:), limdata(1,:), 'Color', limgray, 'LineWidth', 3)
contour(rg,zg,psizr,[psibry psibry], 'color', dgray, 'linewidth', 2);

[rx, zx] = isoflux_xpFinder(psizr, 1.3 ,-1.1, rg, zg);
plot(rx, zx, 'X', 'markersize', 18, 'color', 'b', 'linewidth', 5)


[romp, i] = max(rbbbs);
zomp = zbbbs(i);
rlevels = romp + [0.004:0.004:0.03];
zlevels = zomp*ones(size(rlevels));
psilevels = bicubicHermite(rg,zg,psizr, rlevels, zlevels);
contour(rg,zg,psizr,[psilevels psilevels], 'color', lgray, 'linewidth', 0.5);



% PERFECT
axes(ax(2))
hold on
load('eq_perf')
struct_to_ws(eq.gdata);
[psizr, rg, zg] = regrid(rg, zg, psizr, 300, 300);

plot(limdata(2,:), limdata(1,:), 'Color', limgray, 'LineWidth', 3)
contour(rg,zg,psizr,[psibry psibry], 'color', dgray, 'linewidth', 2);

rx = 1.17;
zx = -1.25;
plot(rx, zx, 'X', 'markersize', 18, 'color', 'b', 'linewidth', 5);

plot(rx, zx, 'X', 'color', 'b', 'linewidth', 4)

[romp, i] = max(rbbbs);
zomp = zbbbs(i);
rlevels = romp + [0.004:0.004:0.03];
zlevels = zomp*ones(size(rlevels));
psilevels = bicubicHermite(rg,zg,psizr, rlevels, zlevels);
contour(rg,zg,psizr,[psilevels psilevels], 'color', lgray, 'linewidth', 0.5);


%%
% SNOW PLUS
axes(ax(3))
cla
hold on
load('eq_sfp')
struct_to_ws(eq.gdata);
[psizr, rg, zg] = regrid(rg, zg, psizr, 300, 300);


plot(limdata(2,:), limdata(1,:), 'Color', limgray, 'LineWidth', 3)
contour(rg,zg,psizr,[psibry psibry], 'color', dgray, 'linewidth', 2);

[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);
plot([rxP rxS], [zxP zxS], 'X', 'markersize', 10, 'color', 'b', 'linewidth', 3);
contour(rg,zg,psizr,[psixS psixS], 'color', dgray, 'linewidth', 2);


[romp, i] = max(rbbbs);
zomp = zbbbs(i);
rlevels = romp + [0.0005:0.0005:0.006];
zlevels = zomp*ones(size(rlevels));
psilevels = bicubicHermite(rg,zg,psizr, rlevels, zlevels);
contour(rg,zg,psizr,[psilevels psilevels], 'color', lgray, 'linewidth', 0.5);



%%
% SNOW MINUS
axes(ax(4))
cla
hold on
load('eq_sfm')
struct_to_ws(eq.gdata);
[psizr, rg, zg] = regrid(rg, zg, psizr, 300, 300);


plot(limdata(2,:), limdata(1,:), 'Color', limgray, 'LineWidth', 3)
contour(rg,zg,psizr,[psibry psibry], 'color', dgray, 'linewidth', 2);

[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);
plot([rxP rxS], [zxP zxS], 'X', 'markersize', 10, 'color', 'b', 'linewidth', 3);

[romp, i] = max(rbbbs);
zomp = zbbbs(i);
rlevels = romp + [0.0005:0.0005:0.006];
zlevels = zomp*ones(size(rlevels));
psilevels = bicubicHermite(rg,zg,psizr, rlevels, zlevels);
contour(rg,zg,psizr,[psilevels psilevels], 'color', lgray, 'linewidth', 0.5);
contour(rg,zg,psizr,[psixS psixS], 'color', dgray, 'linewidth', 2);


%%
if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_snowtypes/fig_snowtypes2.eps';
  saveas(gcf, fn, 'epsc')
end






