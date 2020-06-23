clear
close all

load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155328_sfm/3594/eqs.mat')
% load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp/155328_sfp/2299/eqs.mat')
plot_eq(eqs{end})
struct_to_ws(eqs{end});
[~,~,~,~, psixPL, psixSL] = my_snowfinder(rg, zg, psizr, psibry);
contour(rg,zg,psizr, psixSL - [0:0.01:0.1], 'color', 'g', 'linewidth', 1);
% contour(rg,zg,psizr, psixPL - [0:0.01:0.1], 'color', 'g', 'linewidth', 1);


load('/u/jwai/d3d_snowflake_2020/current/inputs/qperp/qperp_155328.mat')
[~,k] = min(abs(t-3594));
% [~,k] = min(abs(t-2299));


figure(20)
plot(s, qperp(k,:)')
add_lim_colorcode(20)










