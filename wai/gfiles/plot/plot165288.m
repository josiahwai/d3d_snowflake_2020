geteq = 1;

load('/u/jwai/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')
dir = '/u/jwai/d3d_snowflake_2019_wai/gfiles/165286/kefit';
shot = 165286;
times = [2000 2100 2200 2300 2400 2500];


if geteq
    [eq, neq] = read_eq(shot, times/1000, dir);
end


rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
psi_n0 = linspace(0,1,129);
psi_n = linspace(0,1,300);

for ii = 1:neq
    figure(ii)
    clf
    hold on
    subplot(5,1,1)
    ax = gca;
    ax.Position = [0.25 0.55 0.5 0.5];
    plot_d3d_geo(tok_data_struct)
    axis([0.8 2.2 -1.6 -0.8])
    xlabel(['Time: ' num2str(times(ii)) 'ms'])
    psizr = eq.gdata(ii).psizr;
    psibry = eq.gdata(ii).psibry;
    
    
    contour(rg,zg,psizr,[psibry psibry],'b');
    title(strcat(num2str(times(ii)), 'ms'))
    
    
    pprime = interp1(psi_n0, eq.gdata(ii).pprime, psi_n, 'spline');
    ffprim = interp1(psi_n0, eq.gdata(ii).ffprim, psi_n, 'spline');
    subplot(5,1,4)
    plot(psi_n, ffprim)
    title("FF'")
    xlim([0.8 1])
    subplot(5,1,5)
    plot(psi_n, pprime)
    title("P'")
    xlim([0.8 1])
end
hold off



























