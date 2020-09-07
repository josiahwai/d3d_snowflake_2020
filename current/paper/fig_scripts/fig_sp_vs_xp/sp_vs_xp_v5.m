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
load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_155330.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp.mat')
sim_sfp = sim;

load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')
sim_sfp_sp = sim;

shotlist = 155330; 
i_sp = find(ismember(sim_sfp.shots,shotlist));
i_xp = find(ismember(sim_sfp_sp.shots,shotlist));

% find the shots & times that correspond to each other
% for arbitary example k, iuse(k,1) and iuse(k,2) have the same shot/time
iuse = [];
for i = i_sp
  search_time = sim_sfp.times(i);
  search_shot = sim_sfp.shots(i);
  
  for j = i_xp
    if sim_sfp_sp.times(j) == search_time && sim_sfp_sp.shots(j) == search_shot
      iuse = [iuse; i j];
    end
  end
end


i_sp = iuse(:,1);
i_xp = iuse(:,2);

dxp_sp = sim_sfp.dxp(i_sp,:);
dxp_xp = sim_sfp_sp.dxp(i_xp,:);


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
dist_xp = dist_xp(i);

figure
hold on
bar(dist_sp, 'b')
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





