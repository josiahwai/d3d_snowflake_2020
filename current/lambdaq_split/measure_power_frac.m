function [pf_eq, pf_heatflux] = measure_power_frac(eq, snow, s, qperp, lambdaq)

if ~exist('plotit', 'var'), plotit = 0; end

% EQ-ESTIMATED POWER FRACTION
struct_to_ws(eq);
[rzbbbs] = contourc(rg, zg, psizr, [psibry psibry]);
[romp,k] = max(rzbbbs(1,:));
zomp = rzbbbs(2,k);

dpsi = inf;
while abs(dpsi) > .0002
  [psi, psi_r] = bicubicHermite(rg, zg, psizr, romp, zomp);
  dpsi = psi - snow.psixPL;
  romp = romp - .01 * psi_r * dpsi;
end

rsplit = romp;
dpsi = inf;
while abs(dpsi) > .0002
  [psi, psi_r] = bicubicHermite(rg, zg, psizr, rsplit, zomp);
  dpsi = psi - snow.psixSL;
  rsplit = rsplit - .01 * psi_r * dpsi;
end

drsplit = rsplit - romp;

pf_eq = exp(-drsplit/lambdaq);


% MEASURED POWER FRACTION
q = qperp - median(qperp);
if size(q,1) ~= 1, q = q'; end
q(q<0) = 0;

ishadow = find(s > 160.5 & s < 185);
q(ishadow) = 0;



load('d3d_obj_mks_struct_6565.mat');
limdata = tok_data_struct.limdata;
slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
[rlim, zlim] = calcLimDistanceInv(slimtot - s/100, limdata);

r2 = [1.02 1.15];
r4 = [1.15 1.34];

i2 = find( r2(1) < rlim & rlim < r2(2));
i4 = find( r4(1) < rlim & rlim < r4(2));

P2 = 2 * pi * trapz( rlim(i2), rlim(i2) .* q(i2));
P4 = 2 * pi * trapz( rlim(i4), rlim(i4) .* q(i4));

pf_heatflux = P4 / (P2 + P4);

end
