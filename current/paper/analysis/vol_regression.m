[~,iout1] = rmoutliers(sim.djmax);
[~,iout2] = rmoutliers(sim.dV);

iuse = ~(iout1 | iout2);

dV = sim.dV(iuse)';
dj = sim.djmax(iuse)';

w = pinv(dV)*dj;

dj_pred = dV*w;

figure
scatter(dj, dj_pred)
axis equal







