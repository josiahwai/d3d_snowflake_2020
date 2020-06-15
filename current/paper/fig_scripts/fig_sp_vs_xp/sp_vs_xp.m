% Is constraining via strike point consistent with constraining via x-pt?

% ========
% SETTINGS
% ========
saveit = 0;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';

% =================================
% LOAD AND PLOT X-PTS FROM ALL SIMS
% ================================

d = dir([topdir '/**/*xps*']);

dxp_xp = [];
dxp_sp = [];
t_xp = [];
t_sp = [];


for k = 1:length(d)
  
  dirname = d(k).folder;
  t = str2num(dirname(end-3:end));

  load([dirname '/xps.mat'])  % load x-pts
  
  dxp = xps{end} - xps{1};
  
  if contains(dirname, '_sfp_sp')
    dxp_sp = [dxp_sp; dxp];
    t_sp = [t_sp; t];
  elseif contains(dirname, '_sfp')
    if length(xps) > 10
      dxp_xp = [dxp_xp; dxp];
      t_xp = [t_xp; t];
    end
  end

end


i_sp = ismember(t_sp, t_xp);
i_xp = ismember(t_xp, t_sp);

dxp_sp = dxp_sp(i_sp,:);
dxp_xp = dxp_xp(i_xp,:);
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
scatter(dxp_sp(:,1), dxp_sp(:,3), 'r')
scatter(dxp_sp(:,2), dxp_sp(:,4), 'b')

subplot(2,2,2)
hold on
title('XP Constrained')
scatter(dxp_xp(:,1), dxp_xp(:,3), 'r')
scatter(dxp_xp(:,2), dxp_xp(:,4), 'b')

subplot(2,2,3)
hold on
title('Diff')
scatter( dxp_xp(:,1) - dxp_sp(:,1), dxp_xp(:,3) - dxp_sp(:,3), 'r')
scatter( dxp_xp(:,2) - dxp_sp(:,2), dxp_xp(:,4) - dxp_sp(:,4), 'b')

for i = 1:3
  subplot(2,2,i)
  xline(0);
  yline(0);
  box on
  ymax = 0.15;
  axis([-ymax ymax -ymax ymax])
end


if saveit
  fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_dxp/fig_dxp.eps';
  saveas(gcf, fn, 'epsc')
end





