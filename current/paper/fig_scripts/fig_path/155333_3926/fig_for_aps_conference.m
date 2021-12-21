% ========
% SETTINGS
% ========
ccc
shot = 155353;
time_ms = 3926;
saveit = 0;
plot_efit = 0;
simdir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155333_3926/';

% Load stuff
close all; warning('off','all')
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));
load('d3d_obj_mks_struct_129129.mat')

load([simdir 'eqs.mat'])
load([simdir 'xps.mat'])
load([simdir 'sims.mat'])


% colors
blue = [20 108 191]/255;
orange = [198 68 26]/255;

% Define plot axes
figure(10)
clf
spc = .1;
h1 = 0.5;
h2 = 1-h1-3*spc;

ax1 = axes('Position', [0.16 1.03-h1-spc     0.8   h1]); 
ax2 = axes('Position', [0.16 spc            0.8   h2]); 


set(gcf, 'position',  [669 81 409 574])
box(ax1,'on')
box(ax2,'on')

% =================
% PLOT EQUILIBRIA
% =================
axes(ax1)
hold on

% plot limiter and strike points
rlim = tok_data_struct.limdata(2,:); 
zlim = tok_data_struct.limdata(1,:);
plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)


% plot initial eq
% ...............
eq0 = eqs{end};
struct_to_ws(eq0);
clear xlim ylim
snow0 = analyzeSnowflake(eq0);
struct_to_ws(snow0);
[zxP, zxS] = unpack(snow0.zx);
[rxP, rxS] = unpack(snow0.rx);
[psizr, rg, zg] = regrid(rg,zg,psizr,400,400);

% primary contour
crz = contourc(rg,zg,psizr,[psixPL psixPL]); 

% remove a contour line to make figure clearer
k = ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
crz(:,k) = nan;

plot(crz(1,:), crz(2,:), 'Color', 'k', 'linewidth', 1.5);

% secondary contour
crz2 = contourc(rg,zg,psizr,[psixSL psixSL]); 
k = ~inpolygon(crz2(1,:), crz2(2,:), rlim, zlim);
crz2(:,k) = nan;
plot(crz2(1,:), crz2(2,:), 'Color', 'k', 'linewidth', 1.5);

% tertiary contour
crz = contourc(rg,zg,psizr,[psixSL psixSL]-.015); 
k = ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
k = k | crz(1,:) < rxS;
crz(:,k) = nan;
plot(crz(1,:), crz(2,:), 'Color', 'k', 'linewidth', 1);


plot(rx, zx, 'x','Color','k', 'Markersize', 15, 'LineWidth', 4)

% plot formatting
axis equal
axis([0.95 1.45 -1.45 -1])
yticks(-1.4:0.1:-1)
xticks(1:0.1:1.4)
set(gca,'TickDir','in'); 

rsp = [rSPP(1:2) rSPS([2 4])];
zsp = [zSPP(1:2) zSPS([2 4])];
plot(rsp, zsp, 'o', 'color', 'k', 'Markersize', 2, 'LineWidth', 3)
scatter(rsp, zsp, 40, 'filled', 'markerfacecolor', 'k')

t1 = text(rsp(1) -.05, zsp(1),      'SP1');
t2 = text(rsp(2) -.05,  zsp(2)-.008, 'SP2');
t3 = text(rsp(3) -.04,  zsp(3)-.012, 'SP3');
t4 = text(rsp(4) -.015, zsp(4)-.02,  'SP4');
t5 = text(rxP + 0.015, zxP + 0.013, 'XP1');
t6 = text(rxS + 0.02, zxS + 0.01, 'XP2');
set([t1 t2 t3 t4 t5 t6], 'fontweight', 'bold', 'color', 'k');


xlabel('$\mathrm{R [m]}$', 'interpreter', 'latex','fontsize', 12)
ylabel('$\mathrm{Z [m]}$', 'interpreter', 'latex','fontsize', 12)



%%
% ===========
% PLOT HEAT
% ===========
axes(ax2)
cla
hold on

% fit eich profile
struct_to_ws(sims{1});
sir = sir*100;
qirN = qir * sum(qirmaxN)/sum(qirmax);

ef = eich_fitter(sir, qir, eqs{1}, tok_data_struct);

% plot heat flux
plot(sir,qirN,'k','linewidth',1.5)

% find and plot strike points
ssp = [snow0.sSPP(1) 1e6 nan snow0.sSPS(end)];
ef.ssp = double(ef.ssp);
ef.ssp = [ef.ssp(1:2) nan ef.ssp(3)];

for i = [1 2 4]
  xline(ef.ssp(i)*100, '--', 'Color', 'k', 'linewidth', 1);
  
  text(min(ef.ssp(i), ssp(i))*100 - 8, 0.53, ['SP' num2str(i)], ...
    'fontweight', 'bold')
end


% plot formatting
axis([93 182 -0.05 0.6])


yline(0,'-k');
ylabel('$q^{div}_\perp / \sum{q^{div}_{\perp,peaks}}$', ...
  'interpreter', 'Latex','fontsize', 12) 


xlabel('$\mathrm{S [cm]}$', 'interpreter', 'latex','fontsize', 12)


if saveit
  fn = [simdir 'fig_aps.eps'];
  saveas(gcf, fn, 'epsc')
end
