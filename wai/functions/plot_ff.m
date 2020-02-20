% smooth ff'
plotit = 1;
shot = 165288;
time = 2000;
dir = '/u/jwai/d3d_snowflake_2019_wai/gfiles/165288/kefit';


set(0,'defaultlinelinewidth', 1.5);


eq = read_eq(shot,time/1000,dir);

psi_n0 = linspace(0,1,129);
psi_n = linspace(0,1,300);

p =   interp1(psi_n0, eq.gdata.pres,   psi_n, 'spline');
f =   interp1(psi_n0, eq.gdata.fpol,   psi_n, 'spline');
pp =  interp1(psi_n0, eq.gdata.pprime, psi_n, 'spline');
ffp = interp1(psi_n0, eq.gdata.ffprim, psi_n, 'spline');


[pp1] = fd(p, psi_n);
alpha = mean(pp) / mean(pp1);


[fp] = fd(f, psi_n);
ffp1 = f.*fp;
beta = mean(ffp) / mean(ffp1);


k = find(psi_n < 1.1);


if plotit
    figure()
    clf
    
    subplot(2,1,1)
    hold on
    plot(psi_n, pp, 'r')
    plot(psi_n, alpha*pp1, 'b')
    legend('original', 'interpolated')
    title("P'")
    xlim([0.8 1])
    
    subplot(2,1,2)
    hold on
    plot(psi_n(k), ffp(k), 'r');
    plot(psi_n(k), beta*ffp1(k), 'b');
    title("FF'")
    xlim([0.8 1])
end
























