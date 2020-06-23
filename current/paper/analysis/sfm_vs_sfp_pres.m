clear; close all

% ========
% SETTINGS
% ========
use_sfm = 1;
use_sfp = 0;
use_sfp_sp = 1;

saveit = 0;

% ================================
% Linear regression on jpar, dxp
% ================================
load('all_sim_data.mat')
struct_to_ws(sim);

i_sfp_sp = find(i_snowtype == 1);
i_sfp = find(i_snowtype == 2);
i_sfm = find(i_snowtype == 3);

iUse = [];
if use_sfm, iUse = [iUse i_sfm]; end
if use_sfp, iUse = [iUse i_sfp]; end
if use_sfp_sp, iUse = [iUse i_sfp_sp]; end



% FINAL PRESSURE
figure
hold on
subplot(2,1,1)
plot( psin, presf(i_sfm,:), 'r')
ylabel('Pres_{final}')
xlim([0.8 1])
ylim([0 4]*1e4)
title('Snowminus Final')

subplot(2,1,2)
plot( psin, presf(i_sfp_sp,:), 'g')
xlim([0.8 1])
ylim([0 4]*1e4)
title('Snowplus Final')
xlabel('psi')
ylabel('Pres_{final}')



% Initial PRESSURE
figure
hold on
subplot(2,1,1)
plot( psin, pres0(i_sfm,:), 'r')
ylabel('Pres_{Initial}')
xlim([0.8 1])
ylim([0 4]*1e4)
title('Snowminus Initial')

subplot(2,1,2)
plot( psin, pres0(i_sfp_sp,:), 'g')
xlim([0.8 1])
ylim([0 4]*1e4)
title('Snowplus Initial')
xlabel('psi')
ylabel('Pres_{Initial}')


% Delta PRESSURE
figure
subplot(2,1,1)
hold on
plot( psin, dpres(i_sfm,:), 'r')
plot( psin, mean(dpres(i_sfm,:)), 'k', 'linewidth', 2)
ylabel('\Delta P')
xlim([0.8 1])
ylim([-1 1]*1e4)
title('Snowminus \Delta P')

subplot(2,1,2)
hold on
plot( psin, dpres(i_sfp_sp,:), 'g')
plot( psin, mean(dpres(i_sfp_sp,:)), 'k', 'linewidth', 2)
xlim([0.8 1])
ylim([-1 1]*1e4)
title('Snowplus \Delta P')
xlabel('psi')
ylabel('\Delta P')

% i_sfp_sp = [];
% i_sfm = [];

avg_dpres = mean( mean( dpres([i_sfm i_sfp_sp], psin > 0.8)));
abs_avg_dpres = mean( mean( abs( dpres([i_sfm i_sfp_sp], psin > 0.8))));
avg_pres0 = mean( mean( pres0([i_sfm i_sfp_sp], psin > 0.8)));

avg_dpres / avg_pres0
abs_avg_dpres / avg_pres0






















