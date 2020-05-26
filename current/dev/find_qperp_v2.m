
function [s, qdiv_perp, qdiv_par] = find_qperp_v2( ...
  rdiv, zdiv, tau, qentr_par, chi, limdata, eq)

struct_to_ws(eq);

% extend rdiv and zdiv region
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
sdiv = sLimTot - calcLimDistance(rdiv, zdiv, limdata);
ds = mean(diff(sdiv));

n = abs(floor(0.02 / ds)); % extend by 2cm

sdiv_extend = [[sdiv(1) - n*ds: ds: sdiv(1)], sdiv(2:end-1), ...
  [sdiv(end): ds: sdiv(end) + n*ds]];

[rdiv_extend, zdiv_extend] = calcLimDistanceInv(sLimTot - sdiv_extend, limdata);

psi_div_extend = bicubicHermite(rg,zg,psizr,rdiv_extend,zdiv_extend);
psi_div = bicubicHermite(rg,zg,psizr,rdiv,zdiv);

ndiv_extend = length(sdiv_extend);
ndiv = length(sdiv);

  
% initial condition
qentr_par_extend = [zeros(n,1); qentr_par; zeros(n,1)];
qdiv_par = zeros(size(qentr_par_extend));

for ii = 1:ndiv-1
  
  t = tau(ii);
  alpha = sqrt(4*pi*chi*t);
  
  % Define the initial condition
  q0 = qentr_par(ii);
  psi_div_fine = linspace(psi_div(ii), psi_div(ii+1), 100);
  
  % don't solve heat eqn if t==0
  if t ~= 0
    for jj = 1:ndiv_extend
            
      dpsi_div = psi_div_extend(jj) - psi_div_fine;
      
      int = (q0/alpha)*exp(-dpsi_div.^2/(4*chi*t));
      
      qdiv_par(jj) = qdiv_par(jj) - trapz(psi_div_fine, int);
      
    end
  end
end


% sign is flipped, occurs, for example when psibry < 0
if mean(qdiv_par) < 0, qdiv_par = -qdiv_par; end

% Compute the fieldline incidence angles for solution
thetaB = zeros(ndiv_extend,1);

if ~isfield(eq,'bzero')
  bzero = eq.btsurf;
  rzero = eq.rsurf;
end

for k = 1:ndiv_extend
    thetaB(k) = calcBField_IncidenceAngle(rdiv_extend(k), zdiv_extend(k), ...
      psizr, rg, zg, bzero, rzero, limdata);
end

qdiv_perp = qdiv_par.*sin(thetaB);
s = sdiv_extend;


end





















