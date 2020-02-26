function fExp = calcFluxExpansion(eq, rSP, zSP)

% CALCFLUXEXPANSION
%
%   Calculate the flux expansion at the (strike) point located at (rSP,zSP).
%
% USAGE: calcFluxExpansion.m
%
% INPUTS:
%
%   rSP........(strike point) radial coordinate at which to calculate fExp [m]
%
%   zSP........(strike point) vertical coordinate at which to calculate fExp [m]
%
%   psizr.....matrix with dimensions (nz x nr) containing the magnetic flux
%             at nz vertical by nr radial grid points
%
%   rg........array containing the nr radial grid points
%
%   zg........array containing the nz vertical grid points
%
%   psibry....value of the magnetic flux on the plasma boundary [Wb]
%
%   bzero.....applied toroidal field strength at radius rzero [T]
%
%   rzero.....radius at which the toroidal field strength bzero is given [m]
%
% OUTPUTS: 
%
%   fExp......flux expansion at the point (rSP,zSP)
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 02/21/2019
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 02/21/2019
%
%...............................................................................

% Locate the magnetic axis
rg = eq.rg';
zg = eq.zg;
psizr = eq.psizr;
psibry = eq.psibry;
rzero = eq.rzero;
bzero = eq.bzero;
rmaxis0 = eq.rmaxis;
zmaxis0 = eq.zmaxis;

[~, zmaxis, ~] = isoflux_maxisFinder(psizr, rmaxis0, zmaxis0, rg, zg);

% Determine major radius of primary separatrix on the midplane
 
psimid = interp2(rg, zg, psizr, rg, zmaxis);
 
idxP = find(psimid > psibry, 1, 'last' ); 
 
idxPrimary = [idxP-1 idxP idxP+1 idxP+2];

ppPrimary = spline(rg(idxPrimary), psimid(idxPrimary));

cP = ppPrimary.coefs(2,:);

rootsP = roots([cP(1) cP(2) cP(3) cP(4)-psibry]) + rg(idxP);

[~,idxminP] = min(abs(rootsP-rg(idxP)));

rmidPrimary = rootsP(idxminP);

% Compute the poloidal and total fields at the midplane

[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, rmidPrimary, zmaxis);

BrMid = -1/(2*pi*rmidPrimary)*dpsidz;
BzMid =  1/(2*pi*rmidPrimary)*dpsidr;

BpMid = sqrt(BrMid*BrMid + BzMid*BzMid);

BTMid = (bzero*rzero)/rmidPrimary;

BTotMid = sqrt(BpMid*BpMid + BTMid*BTMid);

% Compute the poloidal and total fields at the (strike) point

[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, rSP, zSP);

BrSP = -1/(2*pi*rSP)*dpsidz;
BzSP =  1/(2*pi*rSP)*dpsidr;

BpSP = sqrt(BrSP*BrSP + BzSP*BzSP);

BTSP = (bzero*rzero)/rSP;

BTotSP = sqrt(BpSP*BpSP + BTSP*BTSP);

% Compute the flux expansion at the (strike) point

fExp = (BpMid/BTotMid)/(BpSP/BTotSP);

end
