% Is constraining via strike point consistent with constraining via x-pt?

% ========
% SETTINGS
% ========
clear 
close all
saveit = 0;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';

% =================================
% LOAD AND PLOT X-PTS FROM ALL SIMS
% ================================
load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_jboot_and_xp/sim6_14.mat')
struct_to_ws(sim);

shotlist = 155328:2:155354;
i_sp = find(i_snowtype == 1 & ismember(shots,shotlist));
i_xp = find(i_snowtype == 2 & ismember(shots,shotlist));
% i_sp = find(i_snowtype == 1 & shots == 155354);
% i_xp = find(i_snowtype == 2 & shots == 155354);

% find the shots & times that correspond to each other
% for arbitary example k, iuse(k,1) and iuse(k,2) have the same shot/time
iuse = [];
for i = i_sp
  search_time = times(i);
  search_shot = shots(i);
  
  for j = i_xp
    if times(j) == search_time && shots(j) == search_shot
      iuse = [iuse; i j];
    end
  end
end



i_sp = iuse(:,1);
i_xp = iuse(:,2);

dxp_sp = dxp(i_sp,:);
dxp_xp = dxp(i_xp,:);


[~,k1] = rmoutliers(dxp_sp);
[~,k2] = rmoutliers(dxp_xp);
k = ~(k1 | k2);
dxp_sp = dxp_sp(k,:);
dxp_xp = dxp_xp(k,:);
ddxp = dxp_xp - dxp_sp;


% strike-pt constrained: average sum-of-distances moved by the two x-pts
dist_sp = sqrt( sum( dxp_sp(:,[1 3]).^2, 2)) + sqrt( sum( dxp_sp(:,[2 4]).^2, 2));
disp(['Strike-pt constrained: avg distance moved = ' num2str(mean(dist_sp)) 'cm'])
median(dist_sp);

% x-pt constrained: average sum-of-distances moved by the two x-pts
dist_xp = sqrt( sum( dxp_xp(:,[1 3]).^2, 2)) + sqrt( sum( dxp_xp(:,[2 4]).^2, 2));
disp(['X-pt constrained: avg distance moved = ' num2str(mean(dist_xp)) 'cm'])
median(dist_xp);

% average difference between the strike-pt and x-pt solutions
dx1 = mean( sqrt( sum( ddxp(:,[1 3]).^2, 2)));
dx2 = mean( sqrt( sum( ddxp(:,[2 4]).^2, 2)));



figure

subplot(2,2,1)
hold on
title('SP Constrained')
scatter(dxp_sp(:,1), dxp_sp(:,3), 15, 'filled', 'r')
scatter(dxp_sp(:,2), dxp_sp(:,4),  15, 'filled', 'b')


subplot(2,2,2)
hold on
title('XP Constrained')
scatter(dxp_xp(:,1), dxp_xp(:,3), 15, 'filled', 'r')
scatter(dxp_xp(:,2), dxp_xp(:,4),  15, 'filled', 'b')

subplot(2,2,3)
hold on
title('Diff')
scatter( dxp_xp(:,1) - dxp_sp(:,1), dxp_xp(:,3) - dxp_sp(:,3),  15, 'filled', 'r')
scatter( dxp_xp(:,2) - dxp_sp(:,2), dxp_xp(:,4) - dxp_sp(:,4),  15, 'filled', 'b')

for i = 1:3
  subplot(2,2,i)
  xline(0);
  yline(0);
  box on
  axis([-.05 .05 -.03 .1])
end
subplot(2,2,1)
legend('Primary', 'Secondary', 'location', 'southeast')


if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sp_vs_xp/fig_sp_vs_xp.eps';
  saveas(gcf, fn, 'epsc')
  savefig(gcf, '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sp_vs_xp/fig_sp_vs_xp')
end





