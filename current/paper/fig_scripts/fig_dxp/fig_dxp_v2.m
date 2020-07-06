close all

load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
dxp_sfm = rmoutliers(sim.dxp)*100;


load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')
dxp_sfp = rmoutliers(sim.dxp)*100;


colorblind_cmap

blue = [20 108 191]/255;
orange = [198 68 26]/255;

figure
subplot(1,2,1)
title('Snowflake Minus')
hold on
scatter(dxp_sfm(:,1), dxp_sfm(:,3), 15, 'v',  'markeredgecolor', blue); %, 'markerfacecolor', blue); 
scatter(dxp_sfm(:,2), dxp_sfm(:,4), 15, 'o',  'markeredgecolor', orange); %, 'markerfacecolor', orange);
axis equal
axis([-0.065 0.065 -0.03 0.1]*100)
xline(0, 'k')
yline(0, 'k')
ylabel('\Delta z [cm]')
xlabel('\Delta r [cm]')
legend('xp1', 'xp2')


subplot(1,2,2)
title('Snowflake Plus')
hold on
scatter(dxp_sfp(:,1), dxp_sfp(:,3), 15, 'v',  'markeredgecolor', blue); %, 'markerfacecolor', blue); 
scatter(dxp_sfp(:,2), dxp_sfp(:,4), 15, 'o',  'markeredgecolor', orange); %, 'markerfacecolor', orange);
axis equal
axis([-0.065 0.065 -0.03 0.1]*100)
xline(0, 'k')
yline(0, 'k')
ylabel('\Delta z [cm]')
xlabel('\Delta r [cm]')
legend('xp1', 'xp2')

set(gcf, 'position', [469 382 726 290])

fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_dxp/fig_dxp';
saveas(gcf, [fn '.eps'], 'epsc')






