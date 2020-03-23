% plot convergence
ccc

load('history')

openfig('sfd_geo_335.fig')

x = history.x;
plot(x(:,1),x(:,3),'k','linewidth',2)
plot(x(:,2),x(:,4),'k','linewidth',2)

plot(x(end,1), x(end,3), 'kx', 'Markersize', 12, 'LineWidth', 3)
plot(x(end,2), x(end,4), 'kx', 'Markersize', 12, 'LineWidth', 3)


figure()
plot(history.fval,'linewidth',1.5)
ylabel('J(r_x,z_x)','fontweight','bold')
xlabel('Iteration')
title('X-Point Convergence')
% scatter(x(end,1),x(end,3),200,'kx','filled')
% scatter(x(end,2),x(end,4),200,'kx','filled')