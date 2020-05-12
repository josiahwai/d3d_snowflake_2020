function plot_sim(sim)

figure(27)
% clf
hold on

% plot final heat flux
struct_to_ws(sim);
qir = qir /nansum(qirmax); % renormalize so that sum(pks) = 1
plot(sir,qir,'k','linewidth',1.5)

s = [sI sX sO];
q = [qI; qX; qO];
q = q / nansum(qmax);  % renormalize
plot(s, q, 'r', 'linewidth', 1.5)

axis([0.8 1.8 0 1.05*max(qirmax)/nansum(qirmax)])

set(gcf,'position', [718 -93 535 190])
end