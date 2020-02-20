
simdir = '/u/pvail/d3d_snowflake_2019/invMap/heatsim/';

%...............................................................................
% Gather the heat flux simulation data for shot 157431

times_157431  = 3500:100:5900;

xpData_157431 = zeros(length(times_157431),4);
spData_157431 = zeros(length(times_157431),21);

for ii = 1:length(times_157431)
    
    time = times_157431(ii);
    
    filename = ['HeatFluxSimulation_157431_' int2str(time) '.mat'];
    load([simdir '157431/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_157431(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_157431(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_157431(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
    
    spData_157431(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_157431(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_157431(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_157431(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_157431(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_157431(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_157431(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_157431(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_157431(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
    
    spData_157431(ii,13) = qmax_SP1;
    spData_157431(ii,14) = qmax_SP2;
    spData_157431(ii,15) = qmax_SP3;
    
    spData_157431(ii,16) = sDiv_qmax_SP1;
    spData_157431(ii,17) = sDiv_qmax_SP2;
    spData_157431(ii,18) = sDiv_qmax_SP3;
    
    spData_157431(ii,19) = FWHM_SP1;
    spData_157431(ii,20) = FWHM_SP2;
    spData_157431(ii,21) = FWHM_SP3;
    
    xpData_157431(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_157431(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_157431(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_157431(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 157432

times_157432  = 4000:100:5700;

xpData_157432 = zeros(length(times_157432),4);
spData_157432 = zeros(length(times_157432),21);

for ii = 1:length(times_157432)
    
    time = times_157432(ii);
    
    filename = ['HeatFluxSimulation_157432_' int2str(time) '.mat'];
    load([simdir '157432/' int2str(time) '/' filename])
 
    HeatFluxMapping_ComputeFWHM    
    
    spData_157432(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_157432(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_157432(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_157432(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_157432(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_157432(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_157432(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_157432(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_157432(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_157432(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_157432(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_157432(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_157432(ii,13) = qmax_SP1;
    spData_157432(ii,14) = qmax_SP2;
    spData_157432(ii,15) = qmax_SP3;
    
    spData_157432(ii,16) = sDiv_qmax_SP1;
    spData_157432(ii,17) = sDiv_qmax_SP2;
    spData_157432(ii,18) = sDiv_qmax_SP3;
    
    spData_157432(ii,19) = FWHM_SP1;
    spData_157432(ii,20) = FWHM_SP2;
    spData_157432(ii,21) = FWHM_SP3;
    
    xpData_157432(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_157432(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_157432(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_157432(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158163

times_158163  = 4000:100:5900;

xpData_158163 = zeros(length(times_158163),4);
spData_158163 = zeros(length(times_158163),21);

for ii = 1:length(times_158163)
    
    time = times_158163(ii);
    
    filename = ['HeatFluxSimulation_158163_' int2str(time) '.mat'];
    load([simdir '158163/' int2str(time) '/' filename])

    HeatFluxMapping_ComputeFWHM    
    
    spData_158163(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158163(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158163(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158163(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158163(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158163(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158163(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158163(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158163(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158163(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158163(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158163(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_158163(ii,13) = qmax_SP1;
    spData_158163(ii,14) = qmax_SP2;
    spData_158163(ii,15) = qmax_SP3;
    
    spData_158163(ii,16) = sDiv_qmax_SP1;
    spData_158163(ii,17) = sDiv_qmax_SP2;
    spData_158163(ii,18) = sDiv_qmax_SP3;
    
    spData_158163(ii,19) = FWHM_SP1;
    spData_158163(ii,20) = FWHM_SP2;
    spData_158163(ii,21) = FWHM_SP3;
    
    xpData_158163(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158163(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158163(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158163(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158168
 
times_158168  = 3800:100:5900;
 
xpData_158168 = zeros(length(times_158168),4);
spData_158168 = zeros(length(times_158168),21);

for ii = 1:length(times_158168)
    
    time = times_158168(ii);
    
    filename = ['HeatFluxSimulation_158168_' int2str(time) '.mat'];
    load([simdir '158168/' int2str(time) '/' filename])
 
     HeatFluxMapping_ComputeFWHM   
    
    spData_158168(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158168(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158168(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158168(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158168(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158168(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158168(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158168(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158168(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158168(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158168(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158168(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158168(ii,13) = qmax_SP1;
    spData_158168(ii,14) = qmax_SP2;
    spData_158168(ii,15) = qmax_SP3;
    
    spData_158168(ii,16) = sDiv_qmax_SP1;
    spData_158168(ii,17) = sDiv_qmax_SP2;
    spData_158168(ii,18) = sDiv_qmax_SP3;
    
    spData_158168(ii,19) = FWHM_SP1;
    spData_158168(ii,20) = FWHM_SP2;
    spData_158168(ii,21) = FWHM_SP3;
    
    xpData_158168(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158168(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158168(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158168(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158169
 
times_158169  = 3800:100:5900;
 
xpData_158169 = zeros(length(times_158169),4);
spData_158169 = zeros(length(times_158169),21);

for ii = 1:length(times_158169)
    
    time = times_158169(ii);
    
    filename = ['HeatFluxSimulation_158169_' int2str(time) '.mat'];
    load([simdir '158169/' int2str(time) '/' filename])
 
    HeatFluxMapping_ComputeFWHM
    
    spData_158169(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158169(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158169(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158169(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158169(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158169(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158169(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158169(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158169(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158169(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158169(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158169(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_158169(ii,13) = qmax_SP1;
    spData_158169(ii,14) = qmax_SP2;
    spData_158169(ii,15) = qmax_SP3;
    
    spData_158169(ii,16) = sDiv_qmax_SP1;
    spData_158169(ii,17) = sDiv_qmax_SP2;
    spData_158169(ii,18) = sDiv_qmax_SP3;
    
    spData_158169(ii,19) = FWHM_SP1;
    spData_158169(ii,20) = FWHM_SP2;
    spData_158169(ii,21) = FWHM_SP3;
    
    xpData_158169(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158169(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158169(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158169(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158955

times_158955  = 3800:100:5900;

xpData_158955 = zeros(length(times_158955),4);
spData_158955 = zeros(length(times_158955),21);

for ii = 1:length(times_158955)
    
    time = times_158955(ii);
    
    filename = ['HeatFluxSimulation_158955_' int2str(time) '.mat'];
    load([simdir '158955/' int2str(time) '/' filename])
  
    HeatFluxMapping_ComputeFWHM
    
    spData_158955(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158955(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158955(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158955(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158955(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158955(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158955(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158955(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158955(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158955(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158955(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158955(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158955(ii,13) = qmax_SP1;
    spData_158955(ii,14) = qmax_SP2;
    spData_158955(ii,15) = qmax_SP3;
    
    spData_158955(ii,16) = sDiv_qmax_SP1;
    spData_158955(ii,17) = sDiv_qmax_SP2;
    spData_158955(ii,18) = sDiv_qmax_SP3;
    
    spData_158955(ii,19) = FWHM_SP1;
    spData_158955(ii,20) = FWHM_SP2;
    spData_158955(ii,21) = FWHM_SP3;
    
    xpData_158955(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158955(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158955(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158955(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158956

times_158956  = 3700:100:5900;

xpData_158956 = zeros(length(times_158956),4);
spData_158956 = zeros(length(times_158956),21);

for ii = 1:length(times_158956)
    
    time = times_158956(ii);
    
    filename = ['HeatFluxSimulation_158956_' int2str(time) '.mat'];
    load([simdir '158956/' int2str(time) '/' filename])
 
    HeatFluxMapping_ComputeFWHM
    
    spData_158956(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158956(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158956(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158956(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158956(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158956(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
   
    spData_158956(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158956(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158956(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158956(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158956(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158956(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_158956(ii,13) = qmax_SP1;
    spData_158956(ii,14) = qmax_SP2;
    spData_158956(ii,15) = qmax_SP3;
    
    spData_158956(ii,16) = sDiv_qmax_SP1;
    spData_158956(ii,17) = sDiv_qmax_SP2;
    spData_158956(ii,18) = sDiv_qmax_SP3;
    
    spData_158956(ii,19) = FWHM_SP1;
    spData_158956(ii,20) = FWHM_SP2;
    spData_158956(ii,21) = FWHM_SP3;
    
    xpData_158956(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158956(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158956(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158956(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158957

times_158957  = 3600:100:4300;

xpData_158957 = zeros(length(times_158957),4);
spData_158957 = zeros(length(times_158957),21);

for ii = 1:length(times_158957)
    
    time = times_158957(ii);
    
    filename = ['HeatFluxSimulation_158957_' int2str(time) '.mat'];
    load([simdir '158957/' int2str(time) '/' filename])
 
    HeatFluxMapping_ComputeFWHM
    
    spData_158957(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158957(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158957(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158957(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158957(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158957(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158957(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158957(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158957(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158957(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158957(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158957(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158957(ii,13) = qmax_SP1;
    spData_158957(ii,14) = qmax_SP2;
    spData_158957(ii,15) = qmax_SP3;
    
    spData_158957(ii,16) = sDiv_qmax_SP1;
    spData_158957(ii,17) = sDiv_qmax_SP2;
    spData_158957(ii,18) = sDiv_qmax_SP3;
    
    spData_158957(ii,19) = FWHM_SP1;
    spData_158957(ii,20) = FWHM_SP2;
    spData_158957(ii,21) = FWHM_SP3;
    
    xpData_158957(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158957(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158957(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158957(ii,4)  = HeatFluxSimulation.zxSL;

end

%...............................................................................
% Gather the heat flux simulation data for shot 158959

times_158959  = 3600:100:5100;

xpData_158959 = zeros(length(times_158959),4);
spData_158959 = zeros(length(times_158959),21);

for ii = 1:length(times_158959)
    
    time = times_158959(ii);
    
    filename = ['HeatFluxSimulation_158959_' int2str(time) '.mat'];
    load([simdir '158959/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_158959(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158959(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158959(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158959(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158959(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158959(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158959(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158959(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158959(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158959(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158959(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158959(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158959(ii,13) = qmax_SP1;
    spData_158959(ii,14) = qmax_SP2;
    spData_158959(ii,15) = qmax_SP3;
    
    spData_158959(ii,16) = sDiv_qmax_SP1;
    spData_158959(ii,17) = sDiv_qmax_SP2;
    spData_158959(ii,18) = sDiv_qmax_SP3;
    
    spData_158959(ii,19) = FWHM_SP1;
    spData_158959(ii,20) = FWHM_SP2;
    spData_158959(ii,21) = FWHM_SP3; 
    
    xpData_158959(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158959(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158959(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158959(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 158960

times_158960  = 3500:100:4500;

xpData_158960 = zeros(length(times_158960),4);
spData_158960 = zeros(length(times_158960),21);

for ii = 1:length(times_158960)
    
    time = times_158960(ii);
    
    filename = ['HeatFluxSimulation_158960_' int2str(time) '.mat'];
    load([simdir '158960/' int2str(time) '/' filename])
 
    HeatFluxMapping_ComputeFWHM
    
    spData_158960(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158960(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158960(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158960(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158960(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158960(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158960(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158960(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158960(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158960(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158960(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158960(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158960(ii,13) = qmax_SP1;
    spData_158960(ii,14) = qmax_SP2;
    spData_158960(ii,15) = qmax_SP3;
    
    spData_158960(ii,16) = sDiv_qmax_SP1;
    spData_158960(ii,17) = sDiv_qmax_SP2;
    spData_158960(ii,18) = sDiv_qmax_SP3;
    
    spData_158960(ii,19) = FWHM_SP1;
    spData_158960(ii,20) = FWHM_SP2;
    spData_158960(ii,21) = FWHM_SP3; 
    
    xpData_158960(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158960(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158960(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158960(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 158963

times_158963  = 3600:100:5900;

xpData_158963 = zeros(length(times_158963),4);
spData_158963 = zeros(length(times_158963),21);

for ii = 1:length(times_158963)
    
    time = times_158963(ii);
    
    filename = ['HeatFluxSimulation_158963_' int2str(time) '.mat'];
    load([simdir '158963/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_158963(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158963(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158963(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158963(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158963(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158963(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158963(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158963(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158963(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158963(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158963(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158963(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
    
    spData_158963(ii,13) = qmax_SP1;
    spData_158963(ii,14) = qmax_SP2;
    spData_158963(ii,15) = qmax_SP3;
    
    spData_158963(ii,16) = sDiv_qmax_SP1;
    spData_158963(ii,17) = sDiv_qmax_SP2;
    spData_158963(ii,18) = sDiv_qmax_SP3;
    
    spData_158963(ii,19) = FWHM_SP1;
    spData_158963(ii,20) = FWHM_SP2;
    spData_158963(ii,21) = FWHM_SP3; 
    
    xpData_158963(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158963(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158963(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158963(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 158965

times_158965  = setdiff(3600:100:5900, 5800);

xpData_158965 = zeros(length(times_158965),4);
spData_158965 = zeros(length(times_158965),21);

for ii = 1:length(times_158965)
    
    time = times_158965(ii);
    
    filename = ['HeatFluxSimulation_158965_' int2str(time) '.mat'];
    load([simdir '158965/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_158965(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_158965(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_158965(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_158965(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_158965(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_158965(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_158965(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_158965(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_158965(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_158965(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_158965(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_158965(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_158965(ii,13) = qmax_SP1;
    spData_158965(ii,14) = qmax_SP2;
    spData_158965(ii,15) = qmax_SP3;
    
    spData_158965(ii,16) = sDiv_qmax_SP1;
    spData_158965(ii,17) = sDiv_qmax_SP2;
    spData_158965(ii,18) = sDiv_qmax_SP3;
    
    spData_158965(ii,19) = FWHM_SP1;
    spData_158965(ii,20) = FWHM_SP2;
    spData_158965(ii,21) = FWHM_SP3; 
    
    xpData_158965(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_158965(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_158965(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_158965(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 159008

times_159008  = 3300:100:5900;

xpData_159008 = zeros(length(times_159008),4);
spData_159008 = zeros(length(times_159008),21);

for ii = 1:length(times_159008)
    
    time = times_159008(ii);
    
    filename = ['HeatFluxSimulation_159008_' int2str(time) '.mat'];
    load([simdir '159008/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_159008(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_159008(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_159008(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_159008(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_159008(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_159008(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_159008(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_159008(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_159008(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_159008(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_159008(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_159008(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_159008(ii,13) = qmax_SP1;
    spData_159008(ii,14) = qmax_SP2;
    spData_159008(ii,15) = qmax_SP3;
    
    spData_159008(ii,16) = sDiv_qmax_SP1;
    spData_159008(ii,17) = sDiv_qmax_SP2;
    spData_159008(ii,18) = sDiv_qmax_SP3;
    
    spData_159008(ii,19) = FWHM_SP1;
    spData_159008(ii,20) = FWHM_SP2;
    spData_159008(ii,21) = FWHM_SP3;
    
    xpData_159008(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_159008(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_159008(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_159008(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 165157

times_165157  = 2400:100:3800;

xpData_165157 = zeros(length(times_165157),4);
spData_165157 = zeros(length(times_165157),21);

for ii = 1:length(times_165157)
    
    time = times_165157(ii);
    
    filename = ['HeatFluxSimulation_165157_' int2str(time) '.mat'];
    load([simdir '165157/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_165157(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_165157(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_165157(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_165157(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_165157(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_165157(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_165157(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_165157(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_165157(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_165157(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_165157(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_165157(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_165157(ii,13) = qmax_SP1;
    spData_165157(ii,14) = qmax_SP2;
    spData_165157(ii,15) = qmax_SP3;
    
    spData_165157(ii,16) = sDiv_qmax_SP1;
    spData_165157(ii,17) = sDiv_qmax_SP2;
    spData_165157(ii,18) = sDiv_qmax_SP3;
    
    spData_165157(ii,19) = FWHM_SP1;
    spData_165157(ii,20) = FWHM_SP2;
    spData_165157(ii,21) = FWHM_SP3;
    
    xpData_165157(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_165157(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_165157(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_165157(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 165286

times_165286 = 3200:100:5900;

xpData_165286 = zeros(length(times_165286),4);
spData_165286 = zeros(length(times_165286),21);

for ii = 1:length(times_165286)
    
    time = times_165286(ii);
    
    filename = ['HeatFluxSimulation_165286_' int2str(time) '.mat'];
    load([simdir '165286/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_165286(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_165286(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_165286(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_165286(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_165286(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_165286(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_165286(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_165286(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_165286(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_165286(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_165286(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_165286(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_165286(ii,13) = qmax_SP1;
    spData_165286(ii,14) = qmax_SP2;
    spData_165286(ii,15) = qmax_SP3;
    
    spData_165286(ii,16) = sDiv_qmax_SP1;
    spData_165286(ii,17) = sDiv_qmax_SP2;
    spData_165286(ii,18) = sDiv_qmax_SP3;
    
    spData_165286(ii,19) = FWHM_SP1;
    spData_165286(ii,20) = FWHM_SP2;
    spData_165286(ii,21) = FWHM_SP3;
    
    xpData_165286(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_165286(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_165286(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_165286(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 166147

times_166147 = setdiff(2500:100:5300, [2800 3200 3900 4000 4400 4900]);

xpData_166147 = zeros(length(times_166147),4);
spData_166147 = zeros(length(times_166147),21);

for ii = 1:length(times_166147)
    
    time = times_166147(ii);
    
    filename = ['HeatFluxSimulation_166147_' int2str(time) '.mat'];
    load([simdir '166147/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_166147(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_166147(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_166147(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_166147(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_166147(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_166147(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
    
    spData_166147(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_166147(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_166147(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_166147(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_166147(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_166147(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
 
    spData_166147(ii,13) = qmax_SP1;
    spData_166147(ii,14) = qmax_SP2;
    spData_166147(ii,15) = qmax_SP3;
    
    spData_166147(ii,16) = sDiv_qmax_SP1;
    spData_166147(ii,17) = sDiv_qmax_SP2;
    spData_166147(ii,18) = sDiv_qmax_SP3;
    
    spData_166147(ii,19) = FWHM_SP1;
    spData_166147(ii,20) = FWHM_SP2;
    spData_166147(ii,21) = FWHM_SP3;
    
    xpData_166147(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_166147(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_166147(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_166147(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for shot 175873

times_175873  = [2700 2800 3100:100:4300 4500 4600 4800:100:5900];

xpData_175873 = zeros(length(times_175873),4);
spData_175873 = zeros(length(times_175873),21);

for ii = 1:length(times_175873)
    
    time = times_175873(ii);
    
    filename = ['HeatFluxSimulation_175873_' int2str(time) '.mat'];
    load([simdir '175873/' int2str(time) '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_175873(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_175873(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_175873(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_175873(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_175873(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_175873(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
  
    spData_175873(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_175873(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_175873(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_175873(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_175873(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_175873(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;

    spData_175873(ii,13) = qmax_SP1;
    spData_175873(ii,14) = qmax_SP2;
    spData_175873(ii,15) = qmax_SP3;
    
    spData_175873(ii,16) = sDiv_qmax_SP1;
    spData_175873(ii,17) = sDiv_qmax_SP2;
    spData_175873(ii,18) = sDiv_qmax_SP3;
    
    spData_175873(ii,19) = FWHM_SP1;
    spData_175873(ii,20) = FWHM_SP2;
    spData_175873(ii,21) = FWHM_SP3;
    
    xpData_175873(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_175873(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_175873(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_175873(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Gather the heat flux simulation data for synthetic equilibria

numeqs = [1:35 37:63];

xpData_SYN = zeros(length(numeqs),4);
spData_SYN = zeros(length(numeqs),21);

for ii = 1:length(numeqs)
    
    if ii < 10
        suffixx = ['0' int2str(numeqs(ii))];
    else
        suffixx = int2str(numeqs(ii));
    end
    
    filename = ['HeatFluxSimulation_' suffixx '.mat'];
    load([simdir 'synthetic/eq' suffixx '/' filename])
    
    HeatFluxMapping_ComputeFWHM
    
    spData_SYN(ii,1)  = HeatFluxSimulation.fitEich_SP1.q0;
    spData_SYN(ii,2)  = HeatFluxSimulation.fitEich_SP2.q0;
    spData_SYN(ii,3)  = HeatFluxSimulation.fitEich_SP3.q0;
   
    spData_SYN(ii,4)  = HeatFluxSimulation.fitEich_SP1.S;
    spData_SYN(ii,5)  = HeatFluxSimulation.fitEich_SP2.S;
    spData_SYN(ii,6)  = HeatFluxSimulation.fitEich_SP3.S;
  
    spData_SYN(ii,7)  = HeatFluxSimulation.fitEich_SP1.s0;
    spData_SYN(ii,8)  = HeatFluxSimulation.fitEich_SP2.s0;
    spData_SYN(ii,9)  = HeatFluxSimulation.fitEich_SP3.s0;
    
    spData_SYN(ii,10) = HeatFluxSimulation.fitEich_SP1.fExp;
    spData_SYN(ii,11) = HeatFluxSimulation.fitEich_SP2.fExp;
    spData_SYN(ii,12) = HeatFluxSimulation.fitEich_SP3.fExp;
    
    spData_SYN(ii,13) = qmax_SP1;
    spData_SYN(ii,14) = qmax_SP2;
    spData_SYN(ii,15) = qmax_SP3;
    
    spData_SYN(ii,16) = sDiv_qmax_SP1;
    spData_SYN(ii,17) = sDiv_qmax_SP2;
    spData_SYN(ii,18) = sDiv_qmax_SP3;
    
    spData_SYN(ii,19) = FWHM_SP1;
    spData_SYN(ii,20) = FWHM_SP2;
    spData_SYN(ii,21) = FWHM_SP3;
    
    xpData_SYN(ii,1)  = HeatFluxSimulation.rxPL;
    xpData_SYN(ii,2)  = HeatFluxSimulation.rxSL;
    xpData_SYN(ii,3)  = HeatFluxSimulation.zxPL;
    xpData_SYN(ii,4)  = HeatFluxSimulation.zxSL;
    
end

%...............................................................................
% Aggregate the data

spData_ALL = [spData_157431; spData_157432; spData_158163; spData_158168; ...
              spData_158169; spData_158955; spData_158956; spData_158957; ...
              spData_158959; spData_158960; spData_158963; spData_158965; ...
              spData_159008; spData_165157; spData_165286; spData_166147; ...
              spData_175873; spData_SYN                                   ... 
];

xpData_ALL = [xpData_157431; xpData_157432; xpData_158163; xpData_158168; ...
              xpData_158169; xpData_158955; xpData_158956; xpData_158957; ...
              xpData_158959; xpData_158960; xpData_158963; xpData_158965; ...
              xpData_159008; xpData_165157; xpData_165286; xpData_166147; ...
              xpData_175873; xpData_SYN                                   ...  
];

spData_ALL(:,13) = 0.20*spData_ALL(:,13);
spData_ALL(:,14) = 0.50*spData_ALL(:,14);
spData_ALL(:,15) = 1.00*spData_ALL(:,15);

spData_ALL(:,19) = 3.75*spData_ALL(:,19);
spData_ALL(:,20) = 1.50*spData_ALL(:,20);
spData_ALL(:,21) = 2.50*spData_ALL(:,21);

clear times_157431  times_157432  times_158163  times_158168
clear times_158169  times_158955  times_158956  times_158957
clear times_158959  times_158960  times_158963  times_158965
clear times_159008  times_165157  times_165286  times_166147
clear times_175873

clear spData_157431 spData_157432 spData_158163 spData_158168
clear spData_158169 spData_158955 spData_158956 spData_158957
clear spData_158959 spData_158960 spData_158963 spData_158965
clear spData_159008 spData_165157 spData_165286 spData_166147
clear spData_175873 spData_SYN

clear xpData_157431 xpData_157432 xpData_158163 xpData_158168
clear xpData_158169 xpData_158955 xpData_158956 xpData_158957
clear xpData_158959 xpData_158960 xpData_158963 xpData_158965
clear xpData_159008 xpData_165157 xpData_165286 xpData_166147
clear xpData_175873 spData_SYN

clear simdir filename time ii
clear HeatFluxSimulation

%...............................................................................
% Plot the data

if 0

figure()

subplot(4,3,1)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,1), 'or')
        else % LFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,1), 'ob')
        end
    end
xlabel('s_SP1 [m]', 'Interpreter', 'none')
ylabel('rxP [m]')

subplot(4,3,2)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,1), 'or')
        else % LFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,1), 'ob')
        end
    end
xlabel('s_SP2 [m]', 'Interpreter', 'none')
ylabel('rxP [m]')

subplot(4,3,3)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,1), 'or')
        else % LFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,1), 'ob')
        end
    end
xlabel('s_SP3 [m]', 'Interpreter', 'none')
ylabel('rxP [m]')

subplot(4,3,4)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,2), 'or')
        else % LFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,2), 'ob')
        end
    end
xlabel('s_SP1 [m]', 'Interpreter', 'none')
ylabel('rxS [m]')

subplot(4,3,5)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,2), 'or')
        else % LFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,2), 'ob')
        end
    end
xlabel('s_SP2 [m]', 'Interpreter', 'none')
ylabel('rxS [m]')

subplot(4,3,6)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,2), 'or')
        else % LFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,2), 'ob')
        end
    end
xlabel('s_SP3 [m]', 'Interpreter', 'none')
ylabel('rxS [m]')

subplot(4,3,7)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,3), 'or')
        else % LFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,3), 'ob')
        end
    end
xlabel('s_SP1 [m]', 'Interpreter', 'none')
ylabel('zxP [m]')

subplot(4,3,8)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,3), 'or')
        else % LFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,3), 'ob')
        end
    end
xlabel('s_SP2 [m]', 'Interpreter', 'none')
ylabel('zxP [m]')

subplot(4,3,9)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,3), 'or')
        else % LFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,3), 'ob')
        end
    end
xlabel('s_SP3 [m]', 'Interpreter', 'none')
ylabel('zxP [m]')

subplot(4,3,10)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,4), 'or')
        else % LFS
            plot(spData_ALL(ii,7), xpData_ALL(ii,4), 'ob')
        end
    end
xlabel('s_SP1 [m]', 'Interpreter', 'none')
ylabel('zxS [m]')

subplot(4,3,11)
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,4), 'or')
        else % LFS
            plot(spData_ALL(ii,8), xpData_ALL(ii,4), 'ob')
        end
    end
xlabel('s_SP2 [m]', 'Interpreter', 'none')
ylabel('zxS [m]')

subplot(4,3,12)    
    hold on
    for ii = 1:size(spData_ALL,1)
        if xpData_ALL(ii,1) > xpData_ALL(ii,2) % HFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,4), 'or')
        else % LFS
            plot(spData_ALL(ii,9), xpData_ALL(ii,4), 'ob')
        end
    end
xlabel('s_SP3 [m]', 'Interpreter', 'none')
ylabel('zxS [m]')

end
