% Find the primary and secondary x-pts in a snowflake. 
% This is more robust than snowFinder.m

% Looks for x-pts using Newton's method, with a large number of starting 
% locations. Identifies the best points. 

function [rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry)

% make a grid of initial conditions, look in the lower inner corner only
e = .05;
rmin = min(rg) + e;
rmax = min(rg) + .7*(max(rg) - min(rg));
zmin = min(zg) + e;
zmax = min(zg) + .3*(max(zg) - min(zg)); 

rx0 = linspace(rmin,rmax,20);
zx0 = linspace(zmin,zmax,20);
[rx0,zx0] = meshgrid(rx0,zx0);

% search for x-pts
rx = []; 
zx = [];
for k = 1:numel(rx0)
  try
    [rx(k), zx(k)] = isoflux_xpFinder(psizr,rx0(k),zx0(k),rg,zg);
  catch
    rx(k) = nan;
    zx(k) = nan;
  end
end
rx(isnan(rx)) = [];
zx(isnan(zx)) = [];

% find the best unique x-pts, within range
iuse = find(zx < zmax & zx > zmin & rx > rmin & rx < rmax); 
rxzx = uniquetol([rx(iuse)' zx(iuse)'], .005, 'ByRows',true);

[psix, psix_r, psix_z] = bicubicHermite(rg,zg,psizr,rxzx(:,1), rxzx(:,2));

% only search among pts that are actually x-pts (grad psi == 0)
tol = .001;
igood = find(abs(psix_r.^2 + psix_z.^2) < tol);

[~,k] = mink( abs(psix(igood)-psibry), 2); % take the best 2, if there's > 2
ibest = igood(k);

% % reorder primary, secondary x-pts if necessary
% if abs(psix(k(2))-psibry) < abs(psix(k(1))-psibry), k = flip(k); end

[psixP, psixS] = unpack(psix(ibest));
[rxP, rxS] = unpack(rxzx(ibest,1));
[zxP, zxS] = unpack(rxzx(ibest,2));

end























