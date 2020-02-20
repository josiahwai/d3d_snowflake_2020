
saveit = 0;

%...............................................
% Load the heat flux data for DIII-D shot 165288

datadir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/dataHF/';

datanam_q = 'qperp_165288_4000.txt';
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

figure(22)
plot(s, qperp, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
hold on
plot([50 210], [0 0], '--k')

axis([80 180 -2 50])
xlabel('s [cm]')
ylabel('Heat Flux W/cm^2')
title([int2str(165288) ': ' int2str(4000) ' ms'])

% % EFIT01

% load HeatFluxSimulation_165288_4000_EFIT01.mat
% 
% qdiv_perp_SP1 = 0.80*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP1;
% qdiv_perp_SP2 = 0.80*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP2;
% qdiv_perp_SP3 = 0.80*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP3;
% 
% sDiv_SP1 = 100*HeatFluxSimulation.sDiv_SP1;
% sDiv_SP2 = 100*HeatFluxSimulation.sDiv_SP2;
% sDiv_SP3 = 100*HeatFluxSimulation.sDiv_SP3;
% 
% hold on
% plot(sDiv_SP1, qdiv_perp_SP1, '-', 'color', [0 1 0], 'LineWidth', 2)
% hold on
% plot(sDiv_SP2, qdiv_perp_SP2, '-', 'color', [0 1 0], 'LineWidth', 2)
% hold on
% plot(sDiv_SP3, qdiv_perp_SP3, '-', 'color', [0 1 0], 'LineWidth', 2)


% CAKE Constrained

load HeatFluxSimulation_165288_4000_EQHF.mat

qdiv_perp_SP1 = 1.00*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP1;
qdiv_perp_SP2 = 1.00*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP2;
qdiv_perp_SP3 = 1.00*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP3;

sDiv_SP1 = 100*HeatFluxSimulation.sDiv_SP1;
sDiv_SP2 = 100*HeatFluxSimulation.sDiv_SP2;
sDiv_SP3 = 100*HeatFluxSimulation.sDiv_SP3;

hold on
plot(sDiv_SP1, qdiv_perp_SP1, '-', 'color', [0.91 0.41 0.17], 'LineWidth', 2)
hold on
plot(sDiv_SP2, qdiv_perp_SP2, '-', 'color', [0.91 0.41 0.17], 'LineWidth', 2)
hold on
plot(sDiv_SP3, qdiv_perp_SP3, '-', 'color', [0.91 0.41 0.17], 'LineWidth', 2)



% CAKE unconstrained

% load HeatFluxSimulation_165288_4000_CAKE.mat
% 
% qdiv_perp_SP1 = 0.20*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP1;
% qdiv_perp_SP2 = 0.50*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP2;
% qdiv_perp_SP3 = 1.00*(1e6)*(1e-2)^2*HeatFluxSimulation.qdiv_perp_SP3;
% 
% sDiv_SP1 = 100*HeatFluxSimulation.sDiv_SP1;
% sDiv_SP2 = 100*HeatFluxSimulation.sDiv_SP2;
% sDiv_SP3 = 100*HeatFluxSimulation.sDiv_SP3;
% 
% hold on
% plot(sDiv_SP1, qdiv_perp_SP1, '-', 'color', [0 0 1], 'LineWidth', 2)
% hold on
% plot(sDiv_SP2, qdiv_perp_SP2, '-', 'color', [0 0 1], 'LineWidth', 2)
% hold on
% plot(sDiv_SP3, qdiv_perp_SP3, '-', 'color', [0 0 1], 'LineWidth', 2)



% 
% 
% labels = {'Infrared', 'EFIT01', 'kEFIT + Infrared'};
% co = {[0 0 0], [0 0 1], [0.91 0.41 0.17]};
% location = 'northwest';
% ls = {'-', '-', '-'};
% mylegend(labels,ls,co,location)
% mylegend({'kEFIT + infrared'}, {'-'}, {[0.91 0.41 0.17]}, 'northwest')


if saveit
    saveas(gcf, 'infrared3.eps', 'epsc')
end


