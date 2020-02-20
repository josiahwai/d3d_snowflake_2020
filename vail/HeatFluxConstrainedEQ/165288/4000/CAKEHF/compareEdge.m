
% Compare edge current in EFIT01, CAKE, and CAKEHF

% EFIT
load('/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/EFIT01/eq_165288_4000_EFIT01.mat')
pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,257), 'spline');
ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,257), 'spline');
jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,257), 'spline');
jpar_EFIT01 = jpar_gs;
save('jpar_EFIT01.mat', 'jpar_EFIT01')
psibar_gs = linspace(0,1,257);
psibar = psibar_gs;
save('psibar.mat', 'psibar')

figure()
hold on
plot(psibar_gs, jpar_gs/1e6, 'color', [1 1 1]*0, 'LineWidth', 2)
grid on
xlabel('Normalized Flux \psi_N')
title('Edge Current')


% CAKE
load('/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKE/eq_165288_4000_CAKE.mat')
pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,257), 'spline');
ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,257), 'spline');
jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,257), 'spline');
jpar_CAKE = jpar_gs;
save('jpar_CAKE.mat', 'jpar_CAKE')
psibar_gs = linspace(0,1,257);

hold on
plot(psibar_gs, jpar_gs/1e6, '-r', 'LineWidth', 2)


% CAKE + irtv
load('/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKEHF/eq_165288_4000_CAKEHF.mat')
pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,257), 'spline');
ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,257), 'spline');
jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,257), 'spline');
jpar_CAKEHF = jpar_gs;
save('jpar_CAKEHF.mat', 'jpar_CAKEHF')
psibar_gs = linspace(0,1,257);

hold on
plot(psibar_gs, jpar_gs/1e6, '-', 'color', [153 194 255]/255, 'LineWidth', 2)
xlim([0.85 1.00])


% NEW CAKE + irtv
load('/u/jwai/d3d_snowflake_2019_wai/hf_constrained_eq/eq_cake_constrained_165288_4000.mat')
pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,257), 'spline');
ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,257), 'spline');
jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,257), 'spline');
jpar_NEWCAKEHF = jpar_gs;
save('jpar_NEWCAKEHF.mat', 'jpar_NEWCAKEHF')
psibar_gs = linspace(0,1,257);

hold on
plot(psibar_gs, jpar_gs/1e6, '-', 'color', [0 51 204]/255, 'LineWidth', 2)
% xlim([0.8 1.00])
xlim([0 1])
ylim([0 1.0])



% legend('EFIT01', 'Pats-cake-run', 'New Cake')
legend('EFIT01', 'kinetic EFIT', 'kinetic EFIT + infrared (OLD)', ...
    'kinetic EFIT + infrared (NEW)', 'location', 'northwest')
ylabel('Parallel Current Density J_{||} [MA/m^2]')
