% ========
% SETTINGS
% ========
clear
saveit = 1;

% =======================
% REGRESSION ON JBOOT SFM
% =======================
close all

load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp.mat')

struct_to_ws(sim);

% regression
djmax = djmax * 1e-6; % [A] -> [MA]
if size(djmax,1) == 1, djmax = djmax'; end
w = pinv(dxp) * djmax;  % [MA/cm]

djmax_pred = dxp*w;

explained_variance = 1 - sum( (djmax- djmax_pred).^2) / ...
  sum( djmax.^2);

mad = mean( abs(sim.djmax) ./ sim.j0max);
md = mean(sim.djmax ./ sim.j0max);
disp(['Mean Absolute Deviation: ' num2str(mad)])
disp(['Mean Deviation: ' num2str(md)])
disp(['Explained Variance: ' num2str(explained_variance)])
disp(w)


% ========
% PLOT SFM
% ========
figure
subplot(2,1,1)
hold on
ymax = max([djmax; djmax_pred]);
ymin = min([djmax; djmax_pred]);


scatter(djmax, djmax_pred, 'r')
title('Regression: Snow Minus', 'fontsize', 12)
% xlabel('\Delta j_{true} [MA]', 'fontsize', 12)
ylabel('\Delta j_{predict} [MA]', 'fontsize', 12)

plot( ymax*[-1 1], ymax*[-1 1], '--k', 'linewidth', 0.75)
axis( 1.05*[ymin ymax ymin ymax])

text( 0.25, 0.07, ['Explained variance: ' num2str(floor(explained_variance*100)) '%'],...
  'units', 'normalized')


% =======================
% REGRESSION ON JBOOT SFP
% =======================
load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp.mat')

struct_to_ws(sim);

if size(djmax,1) == 1, djmax = djmax'; end
djmax = djmax * 1e-6; % [A] -> [MA]
k = find(djmax<5);
djmax = djmax(k);
dxp = dxp(k,:);

% regression
w = pinv(dxp) * djmax;  % [MA/cm]
djmax_pred = dxp*w;

explained_variance = 1 - sum( (djmax- djmax_pred).^2) / ...
  sum( djmax.^2);

mad = mean( abs(sim.djmax) ./ sim.j0max);
md = mean(sim.djmax ./ sim.j0max);
disp(['Mean Absolute Deviation: ' num2str(mad)])
disp(['Mean Deviation: ' num2str(md)])
disp(['Explained Variance: ' num2str(explained_variance)])
disp(w)

% ========
% PLOT SFP
% ========
subplot(2,1,2)
hold on
ymax = max([djmax; djmax_pred]);
ymin = min([djmax; djmax_pred]);

scatter(djmax, djmax_pred, 'g')
title('Regression: Snow Plus', 'fontsize', 12)
xlabel('\Delta j_{true} [MA]', 'fontsize', 12)
ylabel('\Delta j_{predict} [MA]', 'fontsize', 12)

plot( ymax*[-1 1], ymax*[-1 1], '--k', 'linewidth', 0.75)
axis( 1.05*[ymin ymax ymin ymax])

text( 0.25, 0.07, ['Explained variance: ' num2str(floor(explained_variance*100)) '%'],...
  'units', 'normalized')

set(gcf,'position',[684 256 332 520])




if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/analysis/';
  fn = [savedir 'fig_regress.eps'];
  saveas(gcf, fn, 'epsc')
end


