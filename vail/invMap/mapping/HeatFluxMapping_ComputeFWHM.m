
sDiv_SP1 = HeatFluxSimulation.sDiv_SP1;
sDiv_SP2 = HeatFluxSimulation.sDiv_SP2;
sDiv_SP3 = HeatFluxSimulation.sDiv_SP3;

qdiv_perp_SP1 = HeatFluxSimulation.qdiv_perp_SP1;
qdiv_perp_SP2 = HeatFluxSimulation.qdiv_perp_SP2;
qdiv_perp_SP3 = HeatFluxSimulation.qdiv_perp_SP3;

% Determine magnitude and location of heat flux peak at each SP

[qmax_SP1, idx_SP1] = max(qdiv_perp_SP1);
[qmax_SP2, idx_SP2] = max(qdiv_perp_SP2);
[qmax_SP3, idx_SP3] = max(qdiv_perp_SP3);

sDiv_qmax_SP1 = sDiv_SP1(idx_SP1);
sDiv_qmax_SP2 = sDiv_SP2(idx_SP2);
sDiv_qmax_SP3 = sDiv_SP3(idx_SP3);

% Determine FWHM at each SP

idx1 = find(qdiv_perp_SP1 >= qmax_SP1/2, 1, 'first');
idx2 = find(qdiv_perp_SP1 >= qmax_SP1/2, 1, 'last');

FWHM_SP1 = abs(sDiv_SP1(idx1) - sDiv_SP1(idx2));

idx1 = find(qdiv_perp_SP2 >= qmax_SP2/2, 1, 'first');
idx2 = find(qdiv_perp_SP2 >= qmax_SP2/2, 1, 'last');

FWHM_SP2 = abs(sDiv_SP2(idx1) - sDiv_SP2(idx2));

idx1 = find(qdiv_perp_SP3 >= qmax_SP3/2, 1, 'first');
idx2 = find(qdiv_perp_SP3 >= qmax_SP3/2, 1, 'last');

FWHM_SP3 = abs(sDiv_SP3(idx1) - sDiv_SP3(idx2));
