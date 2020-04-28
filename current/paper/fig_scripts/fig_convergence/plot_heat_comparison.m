% -----------
% PLOT HEAT
% -----------
ccc
load('efit_sim')
struct_to_ws(sim);

figure(20)
hold on

plot(sir*100,qir*100,'k','linewidth',1.5)
yline(0,'--k')
xlabel('s [cm]')
ylabel('$\mathsf{Heat \; Flux \;\; q^{div}_\perp \; [W/cm^2]}$', ...
  'interpreter', 'Latex') 

axis([90 170 -5 65])



blue = [20 108 191]/255;
orange = [198 68 26]/255;

plot(sI*100, qI*100, 'Color', blue, 'linewidth', 1.5)
plot(sX*100, qX*100, 'Color', blue, 'linewidth', 1.5)
plot(sO*100, qO*100, 'Color', blue, 'linewidth', 1.5)


efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
efit_eq = read_eq(shot, time_ms/1000, efit_dir);
efit_snow = analyzeSnowflake(efit_eq);

















