close all; 
warning('off','all')
shot = 165288;
time_ms = 4200;
saveit = 0;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);


% colors
blue = [20 108 191]/255;
orange = [198 68 26]/255;


% Define plot axes
figure(10)
spc = .07;
h1 = 0.4;
h2 = 0.18;
h3 = 1 - h1 - h2 - spc*4;

ax1 = axes('Position', [0.13 1.02-h1-spc   0.8   h1]); 
ax2 = axes('Position', [0.13 2*spc+h3   0.8   h2]); 
ax3 = axes('Position', [0.13 spc        0.8   h3]); 

set(gcf, 'position', [1453 62 447 811])





% ----------
% PLOT EFIT
% ----------
axes(ax1)
hold on

efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
efit_eq = read_eq(shot, time_ms/1000, efit_dir);
efit_snow = analyzeSnowflake(efit_eq);
struct_to_ws(efit_snow);

contour(rg,zg,psizr,[psixPL psixPL], 'Color', blue, 'linewidth', 1.5);
contour(rg,zg,psizr,[psixSL psixSL], 'Color', blue, 'linewidth', 1.5);
plot(rx, zx, 'x','Color',blue, 'Markersize', 15, 'LineWidth', 4)




% ------------
% PLOT FINAL EQ
% ------------
load('history')

try 
  load('final_eq')
  load('final_sim')
catch  
  xp = history.x(end,:);
  eq = designeq_ml(xp, shot,time_ms);
  sim = heatsim_ml(eq,shot,time_ms);
  J = measure_cost(sim);  
  save('final_eq','eq')
  save('final_sim', 'sim')
end

snow = analyzeSnowflake(eq);
struct_to_ws(snow);

contour(rg,zg,psizr,[psixPL psixPL], 'Color', orange, 'linewidth', 1.5);
contour(rg,zg,psizr,[psixSL psixSL], 'Color', orange, 'linewidth', 1.5);

plot(rx, zx, 'x','Color',orange, 'Markersize', 15, 'LineWidth', 4)

x = history.x;
plot(x(:,1),x(:,3),'k','linewidth',2.5)
plot(x(:,2),x(:,4),'k','linewidth',2.5)


% --------------------------
% Plot limiter and strike pts
% --------------------------
rlim = tok_data_struct.limdata(2,:); 
zlim = tok_data_struct.limdata(1,:);

k = 72:85;
r = [rlim(end)  1.1   0.95   0.95  1.5   1.5 rlim(72:85)]; 
z = [zlim(end) -0.95 -0.95 -1.45 -1.45 -1.25 zlim(72:85)]; 
fill(r,z,'w')


plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
axis equal
axis([0.95 1.45 -1.45 -0.95])
yticks(-1.4:0.1:-1)
xticks(1:0.1:1.4)
set(gca,'TickDir','out'); 


rsp = [rSPP(1:2) rSPS(4)];
zsp = [zSPP(1:2) zSPS(4)];
plot(rsp, zsp, 'xk', 'Markersize', 8, 'LineWidth', 3)

text(rsp(1) -.045, zsp(1),      'SP1', 'fontweight', 'bold')
text(rsp(2) -.04,  zsp(2)-.013, 'SP2', 'fontweight', 'bold')
text(rsp(3) -.015, zsp(3)-.02,  'SP3', 'fontweight', 'bold')

text(1.25, -1.22, '165288: 4200 ms', 'fontsize', 11, 'fontweight', 'bold')

text(0.97, -1.39, 'EFIT01', 'fontsize', 11, 'Color', blue, 'fontweight', 'bold')

text(0.97, -1.42, 'EFIT01 + IRTV constraint', 'fontsize', 11, 'Color', ...
  orange, 'fontweight', 'bold')


xlabel('$\mathrm{R [m]}$', 'interpreter', 'latex','fontsize', 12)
ylabel('$\mathrm{Z [m]}$', 'interpreter', 'latex','fontsize', 12)


text(1.41,-0.965, 'a', 'fontsize', 16, 'fontweight', 'bold')



% -----------
% PLOT HEAT
% -----------
axes(ax2)

load('efit_sim')
struct_to_ws(sim);

% figure(20)
hold on

plot(sir*100,qir*100,'k','linewidth',1.5)
yline(0,'--k')
xlabel('$\mathrm{S [cm]}$', 'interpreter', 'latex','fontsize', 12)
ylabel('$\mathrm{Heat \; Flux \;\; q^{div}_\perp \; [W/cm^2]}$', ...
  'interpreter', 'Latex','fontsize', 12) 

axis([95 170 -5 65])

blue = [20 108 191]/255;
orange = [198 68 26]/255;

% plot efit heat flux
plot(sI*100, qI*100, 'Color', blue, 'linewidth', 1.5)
plot(sX*100, qX*100, 'Color', blue, 'linewidth', 1.5)
plot(sO*100, qO*100, 'Color', blue, 'linewidth', 1.5)

% xline(efit_snow.sSP_heat(1)*100,'--','Color', blue, 'linewidth', 1);
% xline(efit_snow.sSP_heat(2)*100,'--','Color', blue, 'linewidth', 1);
% xline(efit_snow.sSP_heat(4)*100,'--','Color', blue, 'linewidth', 1);


% plot efit+irtv heat flux
load('final_sim')
struct_to_ws(sim);

f = 80; % normalization
plot(sI*100, qI*f, 'Color', orange, 'linewidth', 1.5)
plot(sX*100, qX*f, 'Color', orange, 'linewidth', 1.5)
plot(sO*100, qO*f, 'Color', orange, 'linewidth', 1.5)

% xline(snow.sSP_heat(1)*100,'--','Color', orange, 'linewidth', 1);
% xline(snow.sSP_heat(2)*100,'--','Color', orange, 'linewidth', 1);
% xline(snow.sSP_heat(4)*100,'--','Color', orange, 'linewidth', 1);

text(164,58, 'b', 'fontsize', 16, 'fontweight', 'bold')


text(102.2, 55, 'SP1', 'fontweight', 'bold')
text(130, 55, 'SP2', 'fontweight', 'bold')
text(151.2, 55, 'SP3', 'fontweight', 'bold')


% -----------
% CONVERGENCE
% -----------
axes(ax3)

% figure(30)
hold on
plot(history.fval / max(history.fval),'linewidth',1.5)

ylabel('$\mathrm{J(r_x,z_x)}$','interpreter', 'latex', 'fontsize', 12)

xlabel('$\mathrm{Iteration}$', 'interpreter', 'latex', 'fontsize', 12)
% title('X-Point Convergence')

text(32.5,0.9, 'c', 'fontsize', 17, 'fontweight', 'bold')


box(ax1,'on')
box(ax2,'on')
box(ax3,'on')


if saveit
  saveas(gcf, 'convergence.eps', 'epsc')
end
