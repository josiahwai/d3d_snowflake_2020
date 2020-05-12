close all

% load eqs, sims, xps
N = length(sims);

plot_eq(eqs{1})
axis([1.0 1.5 -1.4 -0.9])

c = flip(cool);
for k = 1:N
  rx = xps{k}(1:2);
  zx = xps{k}(3:4);
  scatter(rx,zx,'filled', 'markerfacecolor', c(floor(k/N*length(c)),:))
end

for k = 1:N
  plotsim(sims{k})
  title(num2str(k))  
%   pause
end