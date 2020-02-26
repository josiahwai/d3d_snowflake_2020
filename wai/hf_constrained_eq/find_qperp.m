% ADD INFO PLEASE

% qdiv_perp: perpendicular heat flux at the strike pt from solving diffusion equation
% sDiv_SP: distance along limiter to heat flux level

function [qdiv_perp, sDiv, qdiv_par] = ...
    find_qperp(nRegion, iRegion, psiperp, tau, q_par_midplane, frad, ...
    chi, psiSOL, rdiv, zdiv, limdata, psizr, rg, zg, bzero, rzero, ...
    sLimTot)


% initialize
psiDivSP = linspace(max(psiperp)+.005, min(psiperp)-.005, nRegion);
qdiv_par = zeros(length(psiDivSP),1);
        
for ii = 1:nRegion-1
    
    t = tau(ii);
    alpha = sqrt(4*pi*chi*t);       
    
    % Define the initial condition
    qentr_par = q_par_midplane(ii)*(1-frad);    
    psiDivFine = linspace(psiSOL(ii), psiSOL(ii+1), 100);
   
   for jj = 1:nRegion % index of points across divertor domain
       
       dpsiDiv = psiDivSP(jj) - psiDivFine;
       
       int = (qentr_par/alpha)*exp(-dpsiDiv.^2/(4*chi*t));
       
       qdiv_par(jj) = qdiv_par(jj) - trapz(abs(psiDivFine), int);
       
   end
end


% Compute the fieldline incidence angles for SP1 solution
thetaB = zeros(nRegion,1);

for k = 1:nRegion
    thetaB(k) = calcBField_IncidenceAngle(rdiv(k), zdiv(k), psizr, rg, zg, ...
        bzero, rzero, limdata);
end

qdiv_perp = qdiv_par.*sin(thetaB);

% distance along limiter to to strike point region
sDiv = sLimTot - calcLimDistance(rdiv, zdiv, limdata);


end