set(0,'DefaultLineLineWidth', 1)
close all
saveit = 1;

load('dxp')
load('dsp')

figure

drsp = dsp(:,2);
iso = [-.0147 .0301] / norm([-.0147 .0301]);
v = dxp(:,[2 4]) * iso';
alpha = v'*drsp / (v'*v);
drsp_pred = ( iso(1) .* dxp(:,2) + iso(2) .* dxp(:,4)) * alpha;


dzsp = dsp(:,3);
iso = [0.0402 -0.0120] / norm( [0.0402 -0.0120]);
v = dxp(:,[2 4]) * iso';
beta = v'*dzsp / (v'*v);
dzsp_pred = ( iso(1) .* dxp(:,2) + iso(2) .* dxp(:,4)) * beta;

% plot z prediction
subplot(1,2,1)
hold on
scatter( dzsp, dzsp_pred, 'r')

plot([-1 1], [-1 1], '--k')
axis( 0.03 * [-1 1 -1 1])
xlabel('\Delta zstrike_{true}')
ylabel('\Delta zstrike_{Predict}')
xticks([-.1:0.02:.1])
yticks([-.1:0.02:.1])

% plot r prediction
subplot(1,2,2)
hold on
scatter( drsp, drsp_pred, 'b')

plot([-1 1], [-1 1], '--k')
axis( 0.07 * [-1 1 -1 1])
xlabel('\Delta rstrike_{true}')
ylabel('\Delta rstrike_{Predict}')
xticks([-.1:0.05:.1])
yticks([-.1:0.05:.1])


set(gcf, 'position', [671 612 577 266])

if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sfp_xp2_correlation/';
  fn = [savedir 'sfp_xp2_correlation.eps'];
  saveas(gcf, fn, 'epsc')
end



















