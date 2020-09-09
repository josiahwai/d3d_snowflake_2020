% ========
% SETTINGS
% ========
clear
saveit = 0;

blue = [20 108 191]/255;
orange = [198 68 26]/255;


% =======================
% REGRESSION ON JBOOT SFM
% =======================
close all

% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm_large_lambdaq.mat')

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
hold on
ymax = max([djmax; djmax_pred]);
ymin = min([djmax; djmax_pred]);


scatter(djmax, djmax_pred, 'markeredgecolor', blue)
title('Edge Current Regression', 'fontsize', 12)
% xlabel('\Delta j_{true} [MA]', 'fontsize', 12)
ylabel('\Delta j_{predict} [MA]', 'fontsize', 12)

text( 0.35, 0.13, ['Snow Minus: R^2=' num2str(floor(explained_variance*100)/100)],...
  'units', 'normalized', 'color', blue, 'fontweight', 'bold')



% =======================
% REGRESSION ON JBOOT SFP
% =======================
load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')

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
scatter(djmax, djmax_pred, 'markeredgecolor', orange)

xlabel('\Delta j_{true} [MA/m^2]', 'fontsize', 12)
ylabel('\Delta j_{predict} [MA/m^2]', 'fontsize', 12)

ymax = 0.7;
ymin = -0.4;
plot( ymax*[-1 1], ymax*[-1 1], '--k', 'linewidth', 0.75)
axis( 1.05*[ymin ymax ymin ymax])


text( 0.35, 0.07, ['Snow Plus:    R^2=' num2str(floor(explained_variance*100)/100)],...
  'units', 'normalized', 'color', orange, 'fontweight', 'bold')

set(gcf,'position',[684 417 400 359])

box on


if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/analysis/';
  fn = [savedir 'fig_regress.eps'];
  saveas(gcf, fn, 'epsc')
end


