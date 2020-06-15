
% Plot settings
% .............
saveit = 0;

close all
figure(27)
set(gcf,'position',[804 223 445 182])
ax1 = gca;
ax1.Position = [0.15 0.25 0.8 0.6];

blue = [20 108 191]/255;
orange = [198 68 26]/255;


% Load and normalize variables
% ............................

% load initial heatflux
struct_to_ws(sims{1});
q0 = [qI; qX; qO];
q0 = q0 / nansum(qmax);  % renormalize
q0(q0<0.001) = nan;
s_qmax0 = s_qmax;
s0 = [sI sX sO];

 
% % load final heatflux
struct_to_ws(sims{end});
qf = [qI; qX; qO];
qf = qf / nansum(qmax);  % renormalize
qf(qf<0.001) = nan;
s_qmaxf = s_qmax;
sf = [sI sX sO];

qir = qir /nansum(qirmax); % renormalize so that sum(pks) = 1

singles2doubles;


% Plot
% ....
axes(ax1)
hold on

plot(s0, q0, 'color', blue, 'linewidth', 1.5)
plot(sf, qf, 'color', orange, 'linewidth', 1.5)
plot(sir,qir,'k','linewidth',1.5)

for i = [1 3]
  try
    text(s_qirmax(i)-.08, 0.49, ['HP' num2str(i)], 'fontweight', 'bold');
    xline(s_qirmax(i),'--k','linewidth', 1);
    xline(s_qmaxf(i),'--','Color', orange, 'linewidth', 1);
    xline(s_qmax0(i),'--','Color', blue, 'linewidth', 1);        
  catch
  end
end



% more plot settings
yline(0,'-k');
set(ax1, 'XTick', 0.8:0.1:2)
set(ax1, 'YTick', 0:0.2:1)
set(ax1, 'box', 'on')
set(ax1, 'XLim', [0.9 1.7])
set(ax1, 'YLim', [-.02 0.65])

xl = xlabel('$\mathrm{S [m]}$', 'interpreter', 'latex','fontsize', 12);
xl.Position(2) = xl.Position(2) - 0.06;

ylabel('$\mathrm{q_\perp / \sum q_{\perp, pks}}$', ...
  'interpreter', 'Latex','fontsize', 12)


if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfm_predict/';
  saveas(gcf, [savedir 'fig_hf.eps'], 'epsc');
end
   



























