function plotsim(sim)

figure(27)
% clf
hold on

% plot final heat flux
struct_to_ws(sim);
qir = qir /nansum(qirmax); % renormalize so that sum(pks) = 1
plot(sir,qir,'k','linewidth',1.5)

s = [sI sX sO];
qI = qI * qirmax(1) * nansum(qmax(2:3)) / qmax(1) / nansum(qirmax(2:3)); 
q = [qI; qX; qO];
q = q / nansum([max(qI) max(qX) max(qO)]);  % renormalize
plot(s, q, 'linewidth', 1.5)

axis([0.8 1.8 0 1.05*max(qmax)/nansum(qmax)])

set(gcf,'position', [718 -93 535 190])
end