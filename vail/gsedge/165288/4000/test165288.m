
% Load the IRTV heat flux data
 
qperp = transpose(importdata('qperp_165288_4000.txt'));
 
s = importdata('s_165288_4000.txt');

% Index data for each SP

idxSP1 = find(s < 115);
idxSP3 = find(s > 145);

idxSP2 = setdiff(1:length(s), [idxSP1; idxSP3]);

% Index data for gap

idxGap = find(s < 170);

% Remove the gap

gap1 = s(idxGap(end));
gap2 = s(idxGap(end)+1);

dgap = gap2 - gap1;

s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;

figure(10)
plot(s, qperp, '-ob', 'LineWidth', 1, 'MarkerSize', 2)

%.................................
% Fit Eich profile to the SP1 data

qEich = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (-(x-s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x-s0))./S) + qBG;

ft = fittype(qEich, 'independent', 'x');

options = fitoptions(ft);
options.StartPoint = [10 1.4 0.2 104 1 10];
options.Lower = [0 0 0 0 0 0];
 
fitEich = fit(s(idxSP1), qperp(idxSP1), ft, options);
 
% Extract the coefficients
 
q0      = fitEich.q0;
S       = fitEich.S;
lambdaQ = fitEich.lambdaQ;
s0      = fitEich.s0;
qBG     = fitEich.qBG;
fExp    = fitEich.fExp;

% Define the optimized Eich profile

s_Eich = linspace(0,300,1000);

qperp_Eich = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (-(s_Eich - s0))./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(s_Eich - s0))./S) + qBG;

hold on
plot(s_Eich, qperp_Eich, '-r', 'LineWidth', 1)

s_SP1 = fitEich.s0;

%.................................
% Fit Eich profile to the SP2 data

qEich = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;

ft = fittype(qEich, 'independent', 'x');

options = fitoptions(ft);
options.StartPoint = [35 1.4 0.2 128 1 10];
options.Lower = [0 0 0 0 0 0];
 
fitEich = fit(s(idxSP2), qperp(idxSP2), ft, options);
 
% Extract the coefficients
 
q0      = fitEich.q0;
S       = fitEich.S;
lambdaQ = fitEich.lambdaQ;
s0      = fitEich.s0;
qBG     = fitEich.qBG;
fExp    = fitEich.fExp;

% Define the optimized Eich profile

s_Eich = linspace(0,300,1000);

qperp_Eich = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (s_Eich - s0)./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (s_Eich - s0)./S) + qBG;

hold on
plot(s_Eich, qperp_Eich, '-r', 'LineWidth', 1)

s_SP2 = fitEich.s0;

%.................................
% Fit Eich profile to the SP3 data

qEich = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;

ft = fittype(qEich, 'independent', 'x');

options = fitoptions(ft);
options.StartPoint = [70 1.4 0.2 158, 1, 10];
options.Lower = [0 0 0 0 0 0 ];
 
fitEich = fit(s(idxSP3), qperp(idxSP3), ft, options);
 
% Extract the coefficients
 
q0      = fitEich.q0;
S       = fitEich.S;
lambdaQ = fitEich.lambdaQ;
s0      = fitEich.s0;
qBG     = fitEich.qBG;
fExp    = fitEich.fExp;

% Define the optimized Eich profile

s_Eich = linspace(0,300,1000);

qperp_Eich = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (s_Eich - s0)./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (s_Eich - s0)./S) + qBG;

hold on
plot(s_Eich, qperp_Eich, '-r', 'LineWidth', 1)

s_SP3 = fitEich.s0;

% % Plot the data and fit
% 
% figure(2)
% plot(s, qperp, '-ob', 'LineWidth', 1, 'MarkerSize', 2)
% 
% hold on
% plot(s_Eich, qperp_Eich, '-g', 'LineWidth', 1)
 
xlim([min(s) max(s)])

xlabel('Distance s [cm]')
ylabel('Heat Flux [W/cm^2]')


% MAPPING

HeatFluxMapping_RegressionTree

rxP_EST = predict(rtree_rxP, 0.01*[s_SP1 s_SP2 s_SP3]);
rxS_EST = predict(rtree_rxS, 0.01*[s_SP1 s_SP2 s_SP3]);
zxP_EST = predict(rtree_zxP, 0.01*[s_SP1 s_SP2 s_SP3]);
zxS_EST = predict(rtree_zxS, 0.01*[s_SP1 s_SP2 s_SP3]);

