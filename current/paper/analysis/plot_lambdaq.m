close all;

shotlist = 155328:2:155354;
shotlist(9) = [];

lamq_sfp = [];
lamq_sfm = [];
k = 0;

for shot = shotlist
  k = k+1;
  
  load(['sol_params2_' num2str(shot) '.mat'])
  struct_to_ws(sol_params);
  
  coeff = coeff(times > 2000, :);
  
  i_sfp = find(coeff(:,5));
  i_sfm = find(~coeff(:,5));
  
  lamq_sfp(k) = mean(coeff(i_sfp,2));
  lamq_sfm(k) = mean(coeff(i_sfm,2));
  
%   lamq_sfp = [lamq_sfp; coeff(i_sfp,2)];
%   lamq_sfm = [lamq_sfm; coeff(i_sfm,2)];

end



subplot(1,2,1)
bar(lamq_sfp)
ylim([0 0.3])
yline(mean(lamq_sfp), 'r')

subplot(1,2,2)
bar(lamq_sfm)
ylim([0 0.3])
yline(nansum(lamq_sfm) / sum(~isnan(lamq_sfm)), 'r');










