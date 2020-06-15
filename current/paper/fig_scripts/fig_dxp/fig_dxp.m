% ========
% SETTINGS
% ========
saveit = 1;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';

% =================================
% LOAD AND PLOT X-PTS FROM ALL SIMS
% ================================

d = dir([topdir '/**/*xps*']);

sfm_dxp = [];
sfp_dxp = [];

for k = 1:length(d)
  
  load([d(k).folder '/xps.mat'])  % load x-pts
  
  if length(xps) >= 3
   dxp = xps{end} - xps{1};
  else
    continue
  end
  
  if contains(d(k).folder, 'sfp_sp')
    sfp_dxp = [sfp_dxp; 100*dxp];
  elseif contains(d(k).folder, 'sfp')
    continue
  else 
    sfm_dxp = [sfm_dxp; 100*dxp];
  end
end

figure

subplot(2,2,1)
hold on
title('Snow-Minus Primary X-Pt')
scatter(sfm_dxp(:,1), sfm_dxp(:,3))
ylabel('\Delta z [cm]')


subplot(2,2,2)
hold on
title('Snow-Minus Secondary X-Pt')
scatter(sfm_dxp(:,2), sfm_dxp(:,4))


subplot(2,2,3)
hold on
title('Snow-Plus Primary X-Pt')
scatter(sfp_dxp(:,1), sfp_dxp(:,3))
xlabel('\Delta r [cm]')
ylabel('\Delta z [cm]')

subplot(2,2,4)
hold on
title('Snow-Plus Secondary X-Pt')
scatter(sfp_dxp(:,2), sfp_dxp(:,4))
xlabel('\Delta r [cm]')




for i = 1:4
  subplot(2,2,i)
  xline(0);
  yline(0);
  
  axis(100*[-0.1 0.1 -0.1 0.1])
end


if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_dxp/fig_dxp';
  saveas(gcf, [fn '.eps'], 'epsc')
  savefig(gcf, fn)
end




    
 