% ========
% SETTINGS
% ========
shot = 155354;
time_ms = 3694;
saveit = 1;
simdir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155354_3694/';
plot_efit = 0;

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
spc = .07;
h1 = 0.4;
h2 = 0.18;
h3 = 1 - h1 - h2 - spc*4;

ax1 = axes('Position', [0.13 1.02-h1-spc   0.8   h1]); 
ax2 = axes('Position', [0.13 2*spc+h3   0.8   h2]); 
ax3 = axes('Position', [0.13 spc        0.8   h3]); 

set(gcf, 'position', [585 -111 453 811])
box(ax1,'on')
box(ax2,'on')
box(ax3,'on')

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
eq0 = eqs{1};
snow0 = analyzeSnowflake(eq0);
struct_to_ws(snow0);
[zxP, zxS] = unpack(snow0.zx);

% primary contour
crz = contourc(rg,zg,psizr,[psixPL psixPL]); 
k = crz(2,:) < zxS;  % remove a contour line to make figure clearer
k = k | ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
crz(:,k) = nan;
plot(crz(1,:), crz(2,:), 'Color', blue, 'linewidth', 1.5);

% secondary contour
crz = contourc(rg,zg,psizr,[psixSL psixSL]); 
k = crz(2,:) > zxP;  
k = k | ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
crz(:,k) = nan;
plot(crz(1,:), crz(2,:), 'Color', blue, 'linewidth', 1.5);

plot(rx, zx, 'x','Color',blue, 'Markersize', 15, 'LineWidth', 4)


% plot final eq
% .............
eq1 = eqs{end};
snow1 = analyzeSnowflake(eq1);
struct_to_ws(snow1);
[zxP, zxS] = unpack(zx);
[rxP, rxS] = unpack(rx);

% primary contour
crz = contourc(rg,zg,psizr,[psixPL psixPL]); 
k = crz(2,:) < zxS;  % remove a contour line to make figure clearer
k = k | ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
crz(:,k) = nan;
plot(crz(1,:), crz(2,:), 'Color', orange, 'linewidth', 1.5);

% secondary contour
crz = contourc(rg,zg,psizr,[psixSL psixSL]); 
k = crz(2,:) > zxP  & crz(1,:) < rxP + .05;
k = k | ~inpolygon(crz(1,:), crz(2,:), rlim, zlim);
crz(:,k) = nan;
plot(crz(1,:), crz(2,:), 'Color', orange, 'linewidth', 1.5);

plot(rx, zx, 'x','Color', orange, 'Markersize', 15, 'LineWidth', 4)


% plot x-pt movement
% ..................
xps{1} = single(xps{1});
xps = reshape(cell2mat(xps),4,[])';

% plot(xps(:,1), xps(:,3),'-ok','linewidth',1.5, 'markersize', 3.5, ...
%   'markerfacecolor', 'k') 
% 
% plot(xps(:,2), xps(:,4),'-ok','linewidth',1.5, 'markersize', 3.5, ...
%   'markerfacecolor', 'k') 

plot(xps(:,1), xps(:,3),'k','linewidth',1.5);
plot(xps(:,2), xps(:,4),'k','linewidth',1.5);


% plot efit x-pts
if plot_efit
  efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
  efit_eq = read_eq(shot, time_ms/1000, efit_dir);
  efit_snow = analyzeSnowflake(efit_eq);
  scatter(efit_snow.rx, efit_snow.zx, 60, 'k', 'filled');
end

% plot formatting
axis equal
axis([0.95 1.45 -1.45 -0.95])
yticks(-1.4:0.1:-1)
xticks(1:0.1:1.4)
set(gca,'TickDir','in'); 

rsp = [rSPP(1:2) rSPS(4)];
zsp = [zSPP(1:2) zSPS(4)];
plot(rsp, zsp, 'xk', 'Markersize', 8, 'LineWidth', 3)

text(rsp(1) -.05, zsp(1),      'SP1', 'fontweight', 'bold')
text(rsp(2) -.05,  zsp(2)-.015, 'SP2', 'fontweight', 'bold')
text(rsp(3) -.015, zsp(3)-.02,  'SP3', 'fontweight', 'bold')

title( [num2str(shot) ': ' num2str(time_ms) 'ms'], 'fontsize', 11, ...
  'fontweight', 'bold')

text(0.97, -1.39, 'EFIT', 'fontsize', 11, 'Color', blue, 'fontweight', 'bold')

text(0.97, -1.42, 'EFIT + IRTV constraint', 'fontsize', 11, 'Color', ...
  orange, 'fontweight', 'bold')


xlabel('$\mathrm{R [m]}$', 'interpreter', 'latex','fontsize', 12)
ylabel('$\mathrm{Z [m]}$', 'interpreter', 'latex','fontsize', 12)


text(1.41,-0.965, 'a', 'fontsize', 16, 'fontweight', 'bold')


% ===========
% PLOT HEAT
% ===========
axes(ax2)
cla
hold on

% fit eich profile
struct_to_ws(sims{1});
sir = sir*100; % [m] to [cm]

ef = eich_fitter(sir, qir, eqs{1}, tok_data_struct);

% plot heat flux
qir = qir /nansum(qirmax); % renormalize so that sum(pks) = 1

s = [flip(sI) sX sO] * 100;
q = [flip(qI); qX; qO];
q = q / nansum(qmax); 

plot(sir,qir,'k','linewidth',1.5)
plot(s, q, 'color', blue, 'linewidth', 1.5)

% find and plot strike points
ssp = [snow0.sSPP(1:2) snow0.sSPS(end)];
ef.ssp = double(ef.ssp);

for i = 1:3
  xline(ef.ssp(i)*100, '--', 'Color', 'k', 'linewidth', 1);
  xline(ssp(i)*100, '--', 'Color', blue, 'linewidth', 1);
  
  text(min(ef.ssp(i), ssp(i))*100 - 6, 0.6, ['SP' num2str(i)], ...
    'fontweight', 'bold')
end


% plot formatting
axis([95 173 -0.05 0.68])

text(0.95, 0.9, 'b', 'units', 'normalized', 'fontsize', 16, ...
  'fontweight', 'bold')

yline(0,'-k');
ylabel('$q^{div}_\perp / \sum{q^{div}_{\perp,peaks}}$', ...
  'interpreter', 'Latex','fontsize', 12) 

set(ax2, 'XTickLabels', [])

text(0.97, 0.72, 'EFIT', 'units', 'normalized', 'fontsize', 9, 'Color', ...
  blue, 'fontweight', 'bold', 'horizontalAlignment', 'right')

text(0.97, 0.6, 'IRTV', 'units', 'normalized', 'fontsize', 9, 'Color', ...
  'k', 'fontweight', 'bold', 'horizontalAlignment', 'right')

% ===========
% PLOT HEAT
% ===========
axes(ax3)
cla
hold on

% fit eich profile
struct_to_ws(sims{end});
sir = sir*100; % [m] to [cm]

ef = eich_fitter(sir, qir, eqs{end}, tok_data_struct);

% plot heat flux
qir = qir /nansum(qirmax); % renormalize so that sum(pks) = 1

s = [flip(sI) sX sO] * 100;
q = [flip(qI); qX; qO];
q = q / nansum(qmax) * 1.2; 

plot(sir,qir,'k','linewidth',1.5)
plot(s, q, 'color', orange, 'linewidth', 1.5)

% find and plot strike points
ssp = [snow1.sSPP(1:2) snow1.sSPS(end)];
ef.ssp = double(ef.ssp);

for i = 1:3
  xline(ef.ssp(i)*100, '--', 'Color', 'k', 'linewidth', 1);
  xline(ssp(i)*100, '--', 'Color', orange, 'linewidth', 1);
  
  text(min(ef.ssp(i), ssp(i))*100 - 6, 0.6, ['SP' num2str(i)], ...
    'fontweight', 'bold')
end

% plot formatting
axis([95 173 -0.05 0.68])

text(0.95, 0.9, 'c', 'units', 'normalized', 'fontsize', 16, ...
  'fontweight', 'bold')

yline(0,'-k');
xlabel('$\mathrm{S [cm]}$', 'interpreter', 'latex','fontsize', 12)
ylabel('$q^{div}_\perp / \sum{q^{div}_{\perp,peaks}}$', ...
  'interpreter', 'Latex','fontsize', 12) 

text(0.97, 0.72, 'EFIT+IRTV', 'units', 'normalized', 'fontsize', 9, 'Color', ...
  orange, 'fontweight', 'bold', 'horizontalAlignment', 'right')

text(0.97, 0.6, 'IRTV', 'units', 'normalized', 'fontsize', 9, 'Color', ...
  'k', 'fontweight', 'bold', 'horizontalAlignment', 'right')



if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/fig_path.eps';
  saveas(gcf, fn, 'epsc')
end
