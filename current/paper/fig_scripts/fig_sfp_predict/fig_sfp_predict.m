load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfp_predict/155328_2099/eqs.mat')

saveit = 1;
shot = 155328;


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
contour(rg,zg,psizr,[psixP psixP], 'Color', blue, 'linewidth', 2);
p = plot(rxP, zxP, 'x', 'color', blue, 'markersize', 15, 'linewidth', 4);


% contour secondary
contour(rg,zg,psizr,[psixS psixS], 'Color', lblue, 'linewidth', 2);
p = plot(rxS, zxS, 'x', 'color', lblue, 'markersize', 15, 'linewidth', 4);


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
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfp_predict/';
  saveas(gcf, [savedir 'fig_sfp.eps'], 'epsc');
end




















