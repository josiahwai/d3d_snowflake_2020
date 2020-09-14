function [qpar, qperp] = heat_diffusion(q0, tau, D, sdiv, rdiv, zdiv, eq, limdata)

if size(q0,2) ~= 1, q0 = q0'; end
if size(sdiv,2) ~= 1, sdiv = sdiv'; end
n = length(sdiv);

% =======================
% SOLVE THE HEAT EQUATION
% =======================
qpar = zeros(n,1);

for i = 1:n
  si = sdiv(i);
  t = tau(i);
  
  % green's function for heat equation
  greens =  1 / sqrt(4*pi*t) * exp( - (sdiv - si).^2 / (4 * D * t));         
  
  qpar(i) = trapz(sdiv, q0 .* greens);
  
end

% ==========================
% CONVERT FROM QPAR TO QPERP
% ==========================
if nanmedian(qpar) < 0, qpar = -qpar; end

if isfield(eq, 'bzero')
  bzero = eq.bzero;
  rzero = eq.rzero;
else
  bzero = eq.btsurf;  
  rzero = eq.rsurf;
end


thetaB = zeros(n,1);
for i = 1:n
    thetaB(i) = calcBField_IncidenceAngle(rdiv(i), zdiv(i), eq.psizr, ...
      eq.rg, eq.zg, eq.bzero, eq.rzero, limdata);        
end

qperp = qpar.*sin(thetaB);   

end






















