
close all

load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
struct_to_ws(sim);
djmax = djmax';
dzmaxis = dV';

[~,iout1] = rmoutliers(djmax);
[~,iout2] = rmoutliers(dzmaxis);

iuse = ~(iout1 | iout2);

dzmaxis = dzmaxis(iuse) - mean(dzmaxis(iuse));
djmax = djmax(iuse) - mean(djmax(iuse));

w = pinv(dzmaxis)*djmax;

dj_pred = dzmaxis*w;

figure
scatter(djmax, dj_pred)
axis equal
axis([-2 2 -2 2]*1e5)



dz = linspace(min(dzmaxis), max(dzmaxis), 100);
dj = w*dz;
figure
hold on
scatter(dzmaxis, djmax)
plot(dz,dj)


