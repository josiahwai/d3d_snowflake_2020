% -----------------------
% PLOT CONVERGENCE IN COST
% -----------------------
saveit = 1;

load('history')

figure()
hold on
plot(history.fval / max(history.fval),'linewidth',1.5)

ylabel('$\bf{J(r_x,z_x)}$','interpreter', 'latex')


xlabel('{Iteration}', 'interpeter', 'latex')
title('X-Point Convergence')

text(30,0.9, 'c', 'fontsize', 16, 'fontweight', 'bold')

if saveit
  saveas(gcf, 'convergence.eps', 'epsc')
end
