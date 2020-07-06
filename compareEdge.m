% Compare edge current
load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/3694/eqs.mat')

% load and interpolate
psin = eqs{2}.psibar;
j0 = eqs{2}.jpar*1e-6; % [A] --> [MA]
jf = eqs{3}.jpar*1e-6;     

j0 = interp1(psin, j0, linspace(0,1,200), 'spline');
jf = interp1(psin, jf, linspace(0,1,200), 'spline');
psin = interp1(psin, psin, linspace(0,1,200), 'spline');


% plot it
figure
hold on
plot(psin, j0, 'b', 'LineWidth', 2)
plot(psin, jf, 'color', [0.91 0.41 0.17], 'LineWidth', 2)

grid on
xlabel('Normalized Flux \psi_N')
ylabel('J_{||} [MA/m^2]')

text(0.73, 0.08, '155354: 3700', 'units', 'normalized', 'fontsize', 10)
ylim([-0.2 1])

title('Current Density') 
xlim([0.85 1.00])
box('on')

set(gcf,'position',[855 330 548 308])
legend('kEFIT', 'kEFIT + IRTV', 'location', 'northwest')

