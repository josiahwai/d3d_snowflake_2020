close all

% load eqs, sims, xps
N = 11; % length(sims);

plot_eq(eqs{1})
axis([1.0 1.5 -1.4 -0.9])

c = flip(cool);
for k = 1:N
  rx = xps{k}(1:2);
  zx = xps{k}(3:4);
  scatter(rx,zx,'filled', 'markerfacecolor', c(floor(k/N*length(c)),:))
end
set(gcf,'position', [614 379 441 321])

plotsim(sims{1})
plotsim(sims{min(N,end)})
set(gcf,'position',[617 135 442 143])

% for k = 1:N
%   plotsim(sims{k})
%   title(num2str(k))  
% %   pause
% end

