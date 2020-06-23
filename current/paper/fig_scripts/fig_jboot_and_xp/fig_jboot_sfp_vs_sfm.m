clear; close all
% ========
% SETTINGS
% ========



% =======================
% Edge current sf+ vs sf-
% =======================

load('sim.mat')
struct_to_ws(sim);

shotlist = [155328 155330 155332 155334 155336 155338 155340 155350 155352 155354];

% find indices
i_sfp_sp =[];
i_sfm = [];

for shot = shotlist
  i = find(i_snowtype == 1 & shots == shot);
  j = find(i_snowtype == 3 & shots == shot);
  
  nslices = min( length(i), length(j));  % use an equal # of slices for sfm and sfp
  
  i = datasample(i, nslices); 
  j = datasample(j, nslices); 
  
  i_sfp_sp = [i_sfp_sp sort(i)];
  i_sfm = [i_sfm sort(j)];
end

% ================
% PLOT FINAL JBOOT
% ================
j_sfp = jfmax(i_sfp_sp);
j_sfm = jfmax(i_sfm);

std_dev = sqrt(var(jfmax([i_sfp_sp i_sfm])))

disp(['Mean jboot in snowplus: ' num2str(floor(mean(j_sfp)))])
disp(['Mean jboot in snowminus: ' num2str(floor(mean(j_sfm)))])


ymax = max([j_sfp j_sfm]);

subplot(1,2,1)
hold on
bar(j_sfp)
ylim([0 ymax])
title('J Final SFD+')
xlabel('Slice index')
ylabel('J [A/m^2]')
yline(mean(j_sfp), 'r', 'linewidth', 2);
text( 0.1, 0.9, ['Mean: ' num2str(floor(mean(j_sfp)))], 'units', 'normalized')


subplot(1,2,2)
hold on
bar(j_sfm)
ylim([0 ymax])
title('J Final SFD-')
xlabel('Slice index')
ylabel('J [A/m^2]')
yline(mean(j_sfm), 'r', 'linewidth', 2);
text( 0.1, 0.9, ['Mean: ' num2str(floor(mean(j_sfm)))], 'units', 'normalized')

set(gcf,'position', [308 261 910 424])




% ================
% PLOT INITIAL JBOOT
% ================
j_sfp = j0max(i_sfp_sp);
j_sfm = j0max(i_sfm);

std_dev = sqrt(var(jfmax([i_sfp_sp i_sfm])))

disp(['Mean jboot in snowplus: ' num2str(floor(mean(j_sfp)))])
disp(['Mean jboot in snowminus: ' num2str(floor(mean(j_sfm)))])


ymax = max([j_sfp j_sfm]);
figure
subplot(1,2,1)
hold on
bar(j_sfp)
ylim([0 ymax])
title('J Initial SFD+')
xlabel('Slice index')
ylabel('J [A/m^2]')
yline(mean(j_sfp), 'r', 'linewidth', 2);
text( 0.1, 0.9, ['Mean: ' num2str(floor(mean(j_sfp)))], 'units', 'normalized')


subplot(1,2,2)
hold on
bar(j_sfm)
ylim([0 ymax])
title('J Initial SFD-')
xlabel('Slice index')
ylabel('J [A/m^2]')
yline(mean(j_sfm), 'r', 'linewidth', 2);
text( 0.1, 0.9, ['Mean: ' num2str(floor(mean(j_sfm)))], 'units', 'normalized')

set(gcf,'position', [308 261 910 424])












