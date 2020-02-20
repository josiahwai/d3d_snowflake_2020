dir = [];
shot = 165288;
time = 4000;
set(0,'defaultlinelinewidth', 1.5);

dir{1} = '/u/jwai/d3d_snowflake_2019_wai/gfiles/165288/from-anthony';
dir{2} = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKEHF';
dir{3} = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/EFIT01';

co = {'b','r'};


figure(311)
clf
for i = 1:2
    eq = read_eq(shot,time/1000,dir{i});
    eqs(i) = eq;
    
    pres = eq.gdata.pres;
    fpol = eq.gdata.fpol;
    pprime = eq.gdata.pprime;
    ffprim = eq.gdata.ffprim;


    psiN = linspace(0,1,129);
    
    subplot(2,2,1) 
    hold on
    plot(psiN, pres, co{i})
    scatter(psiN, pres, 10, co{i}, 'filled')
    xlim([0.8 1.05])
    title('P')
    mylegend({'new-cake','Pats-cake-run'}, {'-','-'}, co)
    
    subplot(2,2,2)
    hold on
    plot(psiN, fpol, co{i})
    scatter(psiN, fpol, 10, co{i}, 'filled')
    xlim([0.8 1.05])
    title('F')
    
    subplot(2,2,3) 
    hold on
    plot(psiN, pprime, co{i})
    scatter(psiN, pprime, 10, co{i}, 'filled')
    xlim([0.8 1.05])
    title("P'")
    xlabel('\psi_N')
    
    subplot(2,2,4)
    hold on
    plot(psiN, ffprim, co{i})
    scatter(psiN, ffprim, 10, co{i}, 'filled')
    xlim([0.8 1.05])
    title("FF'")
    xlabel('\psi_N')
   
    
end



















    
