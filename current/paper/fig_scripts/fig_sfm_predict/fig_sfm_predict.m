load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/3860/eqs.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/3860/sims.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/3860/xps.mat')
% 
% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/155354_3694/eqs.mat')

% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/155354_3727/eqs.mat')
% eqs{1} = eqs{1}.gdata;

saveit = 1;
shot = 155354;


close all
root = '/u/jwai/d3d_snowflake_2020/current/';
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);
blue = [20 108 191]/255;
lblue = [143, 173, 217] / 255;

% PLOT EQUILIBRIUM
struct_to_ws(tok_data_struct); 
rlim = limdata(2,:); 
zlim = limdata(1,:);
struct_to_ws(eqs{1});

[psizr, rg, zg] = regrid(rg, zg, psizr, 400, 400);
[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);

plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on

% contour primary
rz = contourc(rg,zg,psizr,[psixP psixP]);
kIn = inpolygon(rz(1,:), rz(2,:), rlim, zlim);
rz(:,~kIn) = nan;
l = plot(rz(1,:), rz(2,:), 'Color', blue, 'linewidth', 2);
l.Color(4) = 1;
p = plot(rxP, zxP, 'x', 'color', blue, 'markersize', 15, 'linewidth', 4);


% contour secondary
rz = contourc(rg,zg,psizr,[psixS psixS]);
kIn = inpolygon(rz(1,:), rz(2,:), rlim, zlim);
rz(:,~kIn) = nan;
rz(:,rz(2,:) > -1.15) = nan;
l = plot(rz(1,:), rz(2,:), 'color', lblue, 'linewidth', 1.5);
% l.Color(4) = 0.3;
p = plot(rxS, zxS, 'x', 'color', lblue, 'markersize', 15, 'linewidth', 4);
% scatter(rxS, zxS, 50, 'rx', 'linewidth', 4, 'MarkerEdgeAlpha', 0.3)

% fig settings
axis equal
axis([0.91 1.48 -1.45 -0.95])
xlabel('R [m]')
ylabel('Z [m]')
xticks(1:0.1:1.4)
yticks(-1.4:0.1:-1)
set(gcf,'position', [698 412 593 460])


text(0.85,0.95,num2str(shot),'fontsize', 10, 'units', 'normalized', ...
  'fontweight', 'bold')






if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/';
  saveas(gcf, [savedir 'fig_sfm.eps'], 'epsc');
end




















