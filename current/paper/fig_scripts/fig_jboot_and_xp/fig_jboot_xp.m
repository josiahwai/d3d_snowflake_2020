% clear; close all
% ========
% SETTINGS
% ========
use_sfm = 1;
use_sfp = 1;
use_sfp_sp = 0;

saveit = 0;

load('sim_data2.mat')
struct_to_ws(sim);



% ================================
% Linear regression on jpar, dxp
% ================================

% i_sfp_sp = find(i_snowtype == 1);
% i_sfp = find(i_snowtype == 2);
% i_sfm = find(i_snowtype == 3);

i_sfp_sp = find(i_snowtype == 1 & shots == 155354);
i_sfp = find(i_snowtype == 2  & shots == 155354);
i_sfm = find(i_snowtype == 3 & shots == 155354);

iUse = [];
if use_sfm, iUse = [iUse i_sfm]; end
if use_sfp, iUse = [iUse i_sfp]; end
if use_sfp_sp, iUse = [iUse i_sfp_sp]; end



% regression
if size(djmax,1) == 1, djmax = djmax'; end
 
w = pinv( dxp(iUse,:)) * djmax(iUse);
djmax_pred = dxp*w;

explained_variance = 1 - sum( (djmax(iUse) - djmax_pred(iUse)).^2) / ...
  sum( djmax(iUse).^2);

disp(['Explained Variance: ' num2str(explained_variance)])


% plot it
fig1 = figure(1);
hold on


if use_sfm
  scatter( djmax(i_sfm), djmax_pred(i_sfm), 'r') 
end

if use_sfp
  scatter( djmax(i_sfp), djmax_pred(i_sfp), 'b')
end

if use_sfp_sp
  scatter( djmax(i_sfp_sp), djmax_pred(i_sfp_sp), 'g');
end

title('Linear Regression: \Delta j_{predict} = w^T\Deltaxp')
xlabel('\Delta j_{true}')
ylabel('\Delta j_{predict}')


ymax = max([djmax(iUse); djmax_pred(iUse)]);
ymin = min([djmax(iUse); djmax_pred(iUse)]);

plot( ymax*[-1 1], ymax*[-1 1], '--k')
axis( 1.05*[ymin ymax ymin ymax])

legend( 'Snowminus', 'Snowplus')
text( 0.1, 0.9, ['Explained variance: ' num2str(floor(explained_variance*100)) '%'], 'units', 'normalized')
 
% dj1 = jf - cake_j0;
% dj = dj - dj1;

% jboot
figure
hold on
subplot(2,1,1)
plot( psin, dj(i_sfm,:), 'r')
ylabel('\Delta j')
xlim([0.8 1])
ylim([-1.5 1]*1e6)
title('Snowminus')

subplot(2,1,2)
plot( psin, dj(i_sfp_sp,:), 'g')
xlim([0.8 1])
ylim([-1.5 1]*1e6)
title('Snowplus')
xlabel('psi')
ylabel('\Delta j')


% jboot
figure
hold on
subplot(2,1,1)
plot( psin, jf(i_sfm,:), 'r')
ylabel('j_{final}')
% xlim([0.8 1])
ylim([-0.3 2]*1e6)
title('Snowminus')

subplot(2,1,2)
plot( psin, jf(i_sfp_sp,:), 'g')
% xlim([0.8 1])
ylim([-0.3 2]*1e6)
title('Snowplus')
xlabel('psi')
ylabel('j_{final}')



% ===================
% Get some statistics
% ===================
mean( abs(sim.djmax) ./ sim.j0max)


djmax = sim.djmax(iUse);
j0max = sim.j0max(iUse);

[~, ioutl_dj] = rmoutliers(djmax, 'median');
[~, ioutl_j0] = rmoutliers(j0max, 'median');

k = ~(ioutl_dj | ioutl_j0);

djmax = sim.djmax(k);
j0max = sim.j0max(k);

mean( abs(djmax) ./ j0max)


if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_jboot_and_xp/';
  fn = [savedir 'lin_regress.eps'];
  saveas(fig1, fn, 'epsc')
end














