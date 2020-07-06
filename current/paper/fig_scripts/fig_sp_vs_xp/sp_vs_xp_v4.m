% Is constraining via strike point consistent with constraining via x-pt?

% ========
% SETTINGS
% ========
clear 
close all
saveit = 1;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';

% =================================
% LOAD AND PLOT X-PTS FROM ALL SIMS
% ================================
load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_jboot_and_xp/sim6_14.mat')
struct_to_ws(sim);

shotlist = 155328:2:155354;
i_sp = find(i_snowtype == 1 & ismember(shots,shotlist));
i_xp = find(i_snowtype == 2 & ismember(shots,shotlist));

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



% Distances moved
dist_sp = sqrt( sum( dxp_sp(:,[1 3]).^2, 2)) + sqrt( sum( dxp_sp(:,[2 4]).^2, 2));
dist_xp = sqrt( sum( dxp_xp(:,[1 3]).^2, 2)) + sqrt( sum( dxp_xp(:,[2 4]).^2, 2));

[dist_sp, i] = sort(dist_sp, 'descend');
dist_xp = sort(dist_xp,'descend');
% dist_xp = dist_sp(i);

figure
hold on
bar(dist_sp)
bar(dist_xp,'r')
xlabel('Index')
ylabel('Distance [cm]')
title('Snowflake Plus X-point Movements')
legend('Strike-pt constrained', 'X-pt constrained')

set(gcf, 'position', [440 420 516 280])

if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sp_vs_xp/fig_sp_vs_xp.eps';
  saveas(gcf, fn, 'epsc')
  savefig(gcf, '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_sp_vs_xp/fig_sp_vs_xp')
end





