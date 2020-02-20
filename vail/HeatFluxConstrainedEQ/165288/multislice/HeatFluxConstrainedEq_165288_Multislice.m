
% Generate equilibria for DIII-D with IRTV-derived constraints on SFD nulls

ploteq = 1;
plotsp = 0;

times = 3800:100:4000;

%........................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

% Compute total length of the limiter [m]

sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% Configure the plots

for ii = 1:length(times)
    
    if ploteq
        
        figure(ii)
        plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
        hold on
        axis equal
        axis([1.0 1.5 -1.4 -0.9])
        xlabel('R [m]')
        ylabel('Z [m]')
        title([int2str(165288) ': ' int2str(times(ii)) ' ms'])
        
    end
    
end

%....................................................................
% Load the EFIT01 equilibria and determine the strike point locations

eq_dir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/EFIT01';

s_SP1_EFIT01 = zeros(length(times),1);
s_SP2_EFIT01 = zeros(length(times),1);
s_SP3_EFIT01 = zeros(length(times),1);

for ii = 1:length(times)
        
    eq = read_eq(165288, times(ii)/1000, eq_dir);

    rg = eq.gdata.rg; 
    zg = eq.gdata.zg;

    rmaxis = eq.gdata.rmaxis;
    zmaxis = eq.gdata.zmaxis;

    bzero = eq.gdata.bzero;
    rzero = eq.gdata.rzero;

    psizr  = eq.gdata.psizr;
    psimag = eq.gdata.psimag;
    psibry = eq.gdata.psibry;

    [psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);

    % Locate the snowflake in the lower divertor

    rExp   =  1.1500;
    zExp   = -1.2500;
    rhoExp =  0.1000;

    [rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
        snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

    [rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
    [rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

    if abs(psixSL - psibry) < abs(psixPL - psibry)
    
        temp   = psixPL;
        psixPL = psixSL;
        psixSL = temp;
    
        temp1 = rxPL;
        temp2 = zxPL;
        rxPL  = rxSL;
        zxPL  = zxSL;
        rxSL  = temp1;
        zxSL  = temp2;
    
    end

    if ploteq
        
        figure(ii)
        contour(rg, zg, psizr, [psixPL psixPL], 'b', 'LineWidth', 2)
        contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)
    
        plot(rxPL, zxPL, 'xb', 'Markersize', 12, 'LineWidth', 3)
        plot(rxSL, zxSL, 'xb', 'Markersize', 12, 'LineWidth', 3)
        
    end

    % Find the strike points in the EFIT
    
    Rlessthan = limdata(2,:) <  2.0;
    Zlessthan = limdata(1,:) < -0.8;
 
    boxL = Rlessthan & Zlessthan;
    limIdxL = find(boxL ~= 0);

    spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
    spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);

    r_SP1_EFIT01 = spRZLP(1,1); 
    r_SP2_EFIT01 = spRZLP(2,1); 
    r_SP3_EFIT01 = spRZLS(4,1);

    z_SP1_EFIT01 = spRZLP(1,2); 
    z_SP2_EFIT01 = spRZLP(2,2); 
    z_SP3_EFIT01 = spRZLS(4,2);

    if ploteq
        
        plot(r_SP1_EFIT01, z_SP1_EFIT01, 'ob', 'LineWidth', 2)
        plot(r_SP2_EFIT01, z_SP2_EFIT01, 'ob', 'LineWidth', 2)
        plot(r_SP3_EFIT01, z_SP3_EFIT01, 'ob', 'LineWidth', 2)
        
    end

    fExp_SP1_EFIT01 = calcFluxExpansion(r_SP1_EFIT01, z_SP1_EFIT01, psizr, ...
        rg, zg, psibry, bzero, rzero);
    fExp_SP2_EFIT01 = calcFluxExpansion(r_SP2_EFIT01, z_SP2_EFIT01, psizr, ...
        rg, zg, psibry, bzero, rzero);
    fExp_SP3_EFIT01 = calcFluxExpansion(r_SP3_EFIT01, z_SP3_EFIT01, psizr, ...
        rg, zg, psibry, bzero, rzero);

    s_SP1_EFIT01(ii) = sLimTot - calcLimDistance(r_SP1_EFIT01, z_SP1_EFIT01, ...
        limdata); 
    s_SP2_EFIT01(ii) = sLimTot - calcLimDistance(r_SP2_EFIT01, z_SP2_EFIT01, ...
        limdata);
    s_SP3_EFIT01(ii) = sLimTot - calcLimDistance(r_SP3_EFIT01, z_SP3_EFIT01, ...
        limdata);

end

%...............................................
% Load the heat flux data for DIII-D shot 165288

datadir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/dataHF/';

qmax_SP1_v = zeros(length(times),1);
qmax_SP2_v = zeros(length(times),1);
qmax_SP3_v = zeros(length(times),1);

sDiv_qmax_SP1_v = zeros(length(times),1);
sDiv_qmax_SP2_v = zeros(length(times),1);
sDiv_qmax_SP3_v = zeros(length(times),1);

FWHM_SP1_v = zeros(length(times),1);
FWHM_SP2_v = zeros(length(times),1);
FWHM_SP3_v = zeros(length(times),1);

s_SP1_IRTV = zeros(length(times),1);
s_SP2_IRTV = zeros(length(times),1);
s_SP3_IRTV = zeros(length(times),1);

for ii = 1:length(times)
    
    datanam_q = ['qperp_165288_' int2str(times(ii)) '.txt'];
    datanam_s = 's_165288.txt';
    
    qperp = transpose(importdata([datadir datanam_q]));
    s = transpose(importdata([datadir datanam_s]));
    
    % Remove the gap
    
    idxGap = find(s < 170);
    
    % Remove the gap
    
    gap1 = s(idxGap(end));
    gap2 = s(idxGap(end)+1);
    
    dgap = gap2 - gap1;
    
    s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;
    
    % Index data for each SP
    
    idx_SP1 = find(s < 115);
    idx_SP3 = find(s > 145);
    
    idx_SP2 = setdiff(1:length(s), [idx_SP1 idx_SP3]);
    
    qperp_SP1 = qperp(idx_SP1);
    qperp_SP2 = qperp(idx_SP2);
    qperp_SP3 = qperp(idx_SP3);
    
    s_SP1 = s(idx_SP1)';
    s_SP2 = s(idx_SP2)';
    s_SP3 = s(idx_SP3)';
     
    % Determine magnitude and location of heat flux peak at each SP
 
    [qmax_SP1_temp, idx_SP1] = max(qperp_SP1);
    [qmax_SP2_temp, idx_SP2] = max(qperp_SP2);
    [qmax_SP3_temp, idx_SP3] = max(qperp_SP3);
    
    qmax_SP1_v(ii) = qmax_SP1_temp*(1e-6)*(100)^2;
    qmax_SP2_v(ii) = qmax_SP2_temp*(1e-6)*(100)^2;
    qmax_SP3_v(ii) = qmax_SP3_temp*(1e-6)*(100)^2;
    
    sDiv_qmax_SP1_v(ii) = 0.01*s_SP1(idx_SP1);
    sDiv_qmax_SP2_v(ii) = 0.01*s_SP2(idx_SP2);
    sDiv_qmax_SP3_v(ii) = 0.01*s_SP3(idx_SP3);
 
    % Determine FWHM at each SP
 
    idx11 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'first');
    idx21 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'last'); 
    
    FWHM_SP1_v(ii) = 0.01*abs(s_SP1(idx11) - s_SP1(idx21));
 
    idx12 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'first');
    idx22 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'last'); 
    
    FWHM_SP2_v(ii) = 0.01*abs(s_SP2(idx12) - s_SP2(idx22));
 
    idx13 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'first');
    idx23 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'last');
 
    FWHM_SP3_v(ii) = 0.01*abs(s_SP3(idx13) - s_SP3(idx23));
 
    %....................................................
    % Fit Eich profile to IRTV data for each strike point
 
    qEich_Inner = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
        (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (-(x-s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x-s0))./S) + qBG;
    
    qEich_Outer = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
        (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;
    
    qGauss = @(q0, s0, sigma, qBG, x) q0*exp(-(x-s0).^2/sigma) + qBG;
    
    ftInner = fittype(qEich_Inner, 'problem', 'fExp', 'independent', 'x');
    ftOuter = fittype(qEich_Outer, 'problem', 'fExp', 'independent', 'x');
    ftGauss = fittype(qGauss,                         'independent', 'x');
    
    optionsInner = fitoptions(ftInner);
    optionsInner.Lower = [0 0 0 0 0];
    
    optionsOuter = fitoptions(ftOuter);
    optionsOuter.Lower = [0 0 0 0 0];
    
    optionsGauss = fitoptions(ftGauss);
    optionsGauss.Lower = [0 0 0 0];
    
    %........................
    % Eich profile fit to SP1
    
    [maxq, idx_maxq] = max(qperp_SP1);
    
    optionsInner.StartPoint = [2*maxq 1.4 0.2 s_SP1(idx_maxq) 0];
    
    fitEich_SP1 = fit(s_SP1, qperp_SP1, ftInner, 'problem', fExp_SP1_EFIT01, ...
        optionsInner);
        
    %.........................
    % Gauss profile fit to SP2
    
    [maxq, idx_maxq] = max(qperp_SP2);
    
    optionsGauss.StartPoint = [2*maxq s_SP2(idx_maxq) 1 0]; 
    
    fitGauss_SP2 = fit(s_SP2, qperp_SP2, ftGauss, optionsGauss);
    
    %........................
    % Eich profile fit to SP3
    
    [maxq, idx_maxq] = max(qperp_SP3);
    optionsOuter.StartPoint = [2*maxq 1.4 0.2 s_SP3(idx_maxq) 0];
    
    fitEich_SP3 = fit(s_SP3, qperp_SP3, ftOuter, 'problem', fExp_SP3_EFIT01, ...
        optionsOuter);
        
    s_SP1_IRTV(ii) = 0.01*fitEich_SP1.s0;
    s_SP2_IRTV(ii) = 0.01*fitGauss_SP2.s0;
    s_SP3_IRTV(ii) = 0.01*fitEich_SP3.s0;
    
    if ploteq
        
        sInv1 = sLimTot - s_SP1_IRTV(ii);
        sInv2 = sLimTot - s_SP2_IRTV(ii);
        sInv3 = sLimTot - s_SP3_IRTV(ii);
        
        [r_SP1_IRTV, z_SP1_IRTV] = calcLimDistanceInv(sInv1, limdata);
        [r_SP2_IRTV, z_SP2_IRTV] = calcLimDistanceInv(sInv2, limdata);
        [r_SP3_IRTV, z_SP3_IRTV] = calcLimDistanceInv(sInv3, limdata);
        
        figure(ii)
        plot(r_SP1_IRTV, z_SP1_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2)
        plot(r_SP2_IRTV, z_SP2_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2)
        plot(r_SP3_IRTV, z_SP3_IRTV, 'xk', 'MarkerSize', 8, 'LineWidth', 2) 
        
    end

end

%....................................................
% Predict the null locations using the heat flux data

rxP_PRED = zeros(length(times),1);
rxS_PRED = zeros(length(times),1);
zxP_PRED = zeros(length(times),1);
zxS_PRED = zeros(length(times),1);

% HeatFluxMapping_RegressionTree_v6
HeatFluxMapping_RegressionTree_v1

for ii = 1:length(times)
          
%     rxP_PRED(ii) = predict(rtree_rxP, [                             ...
%         qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
%         sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
%         FWHM_SP1_v(ii)      FWHM_SP2_v(ii)      FWHM_SP3_v(ii)      ...
%     ]);
% 
%     rxS_PRED(ii) = predict(rtree_rxS, [                             ...
%         qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
%         sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
%         FWHM_SP1_v(ii)      FWHM_SP2_v(ii)      FWHM_SP3_v(ii)      ...
%     ]);
% 
%     zxP_PRED(ii) = predict(rtree_zxP, [                             ...
%         qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
%         sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
%         FWHM_SP1_v(ii)      FWHM_SP2_v(ii)      FWHM_SP3_v(ii)      ...
%     ]);
% 
%     zxS_PRED(ii) = predict(rtree_zxS, [                             ...
%         qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
%         sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
%         FWHM_SP1_v(ii)      FWHM_SP2_v(ii)      FWHM_SP3_v(ii)      ...
%     ]);

    rxP_PRED(ii) = predict(rtree_rxP, [                             ...
        sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
    ]);

    rxS_PRED(ii) = predict(rtree_rxS, [                             ...
        sDiv_qmax_SP1_v(ii) sDiv_qmax_SP2_v(ii) sDiv_qmax_SP3_v(ii) ...
    ]);

    zxP_PRED(ii) = predict(rtree_zxP, [                             ...
        qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
    ]);

    zxS_PRED(ii) = predict(rtree_zxS, [                             ...
        qmax_SP1_v(ii) qmax_SP2_v(ii) qmax_SP3_v(ii) ...
    ]);

    if ploteq

        figure(ii)
        plot(rxP_PRED(ii), zxP_PRED(ii), 'x', 'color', [0.91 0.41 0.17], ...
            'Markersize', 12, 'LineWidth', 3)
        hold on
        plot(rxS_PRED(ii), zxS_PRED(ii), 'x', 'color', [0.91 0.41 0.17], ...
            'Markersize', 12, 'LineWidth', 3)
        
    end
    
end

%..................................................................
% Design a Grad-Shafranov equilibrium with constraints on the nulls
 
addpath(genpath('/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/gatools'))

s_SP1_EQHF = zeros(length(times),1);
s_SP2_EQHF = zeros(length(times),1);
s_SP3_EQHF = zeros(length(times),1);

for ii = 1:length(times)
    
    gsdesign_d3d_165288_EQHF
    
    rg = eq.rg; 
    zg = eq.zg;
    
    psizr  = eq.psizr;
    psibry = eq.psibry;
    
    [psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);
    
    % Locate the snowflake in the lower divertor
    
    rExp   =  1.1500;
    zExp   = -1.2500;
    rhoExp =  0.1000;
    
    [rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
        snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);
    
    [rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
    [rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);
    
    if abs(psixSL - psibry) < abs(psixPL - psibry)
        
        temp   = psixPL;
        psixPL = psixSL;
        psixSL = temp;
        
        temp1 = rxPL;
        temp2 = zxPL;
        rxPL  = rxSL;
        zxPL  = zxSL;
        rxSL  = temp1;
        zxSL  = temp2;
    
    end
    
    if ploteq
        
        figure(ii)
        contour(rg, zg, psizr, [psixPL psixPL], 'color', [0.91 0.41 0.17], ...
            'LineWidth', 2)
        contour(rg, zg, psizr, [psixSL psixSL], 'color', [0.91 0.41 0.17], ...
            'LineWidth', 1)
        
        plot(rxPL, zxPL, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
            'LineWidth', 3)
        plot(rxSL, zxSL, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
            'LineWidth', 3)
        
    end
    
    % Find the strike points in the kinetic EFIT
    
    Rlessthan = limdata(2,:) <  2.0;
    Zlessthan = limdata(1,:) < -0.8;
    
    boxL = Rlessthan & Zlessthan;
    limIdxL = find(boxL ~= 0);
    
    spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
    spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);
    
    r_SP1 = spRZLP(1,1); r_SP2 = spRZLP(2,1); r_SP3 = spRZLS(4,1);
    z_SP1 = spRZLP(1,2); z_SP2 = spRZLP(2,2); z_SP3 = spRZLS(4,2);
    
    % Compute distance to strike points [m]
    
    s_SP1_EQHF(ii) = sLimTot - calcLimDistance(r_SP1, z_SP1, limdata); 
    s_SP2_EQHF(ii) = sLimTot - calcLimDistance(r_SP2, z_SP2, limdata);
    s_SP3_EQHF(ii) = sLimTot - calcLimDistance(r_SP3, z_SP3, limdata);

end

figure(ii)
hold on

plot(times, 1000*s_SP1_IRTV, '-ok', 'LineWidth', 2)
plot(times, 1000*s_SP2_IRTV, '-ok', 'LineWidth', 2)
plot(times, 1000*s_SP3_IRTV, '-ok', 'LineWidth', 2)

plot(times, 1000*s_SP1_EFIT01, '-ob', 'LineWidth', 2)
plot(times, 1000*s_SP2_EFIT01, '-ob', 'LineWidth', 2)
plot(times, 1000*s_SP3_EFIT01, '-ob', 'LineWidth', 2)

plot(times, 1000*s_SP1_EQHF, '-o', 'color', [0.91 0.41 0.17], 'LineWidth', 2)
plot(times, 1000*s_SP2_EQHF, '-o', 'color', [0.91 0.41 0.17], 'LineWidth', 2)
plot(times, 1000*s_SP3_EQHF, '-o', 'color', [0.91 0.41 0.17], 'LineWidth', 2)

xlabel('Time [ms]')
ylabel('Strike Point Location [cm]')
title('165288')
