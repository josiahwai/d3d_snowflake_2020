% Compare edge current
load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp/155354_sfp/3561/eqs.mat')

% load and interpolate
psin = eqs{2}.psibar;
j0 = eqs{2}.jpar*1e-6; % [A] --> [MA]
jf = eqs{3}.jpar*1e-6;     

j0 = interp1(psin, j0, linspace(0,1,200), 'spline');
jf = interp1(psin, jf, linspace(0,1,200), 'spline');
psin = interp1(psin, psin, linspace(0,1,200), 'spline');

blue = [20 108 191]/255;
orange = [198 68 26]/255;

% plot it
figure
hold on
plot(psin, j0, 'color', blue, 'LineWidth', 2)
plot(psin, jf, 'color', orange, 'LineWidth', 2)

grid on
set(gca,'gridAlpha',0.2)
xlabel('Normalized Flux \psi_N')
ylabel('J_{||} [MA/m^2]')

text(0.73, 0.08, '155354: 3560', 'units', 'normalized', 'fontsize', 10)
ylim([-0.2 1])

title('Current Density') 
xlim([0.85 1.00])
box('on')

set(gcf,'position',[855 330 548 308])
legend('kEFIT', 'kEFIT + IRTV', 'location', 'northwest')

fn = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_dj_155354_3700/fig_dj.eps';
saveas(gcf,fn,'epsc')
