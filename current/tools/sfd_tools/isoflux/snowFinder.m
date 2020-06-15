function [rxP, zxP, rxS, zxS, rSnow, zSnow, dR, dZ, rho, theta] = ...
    snowFinder(psizr, rExp, zExp, rhoExp, rg, zg, psibry)
%
% SNOWFINDER
%
%   Locate the two magnetic nulls of a snowflake configuration by using a
%   third-order expansion of the magnetic flux function to compute the
%   field structure in the divertor region.
%
%   Described in D.D. Ryutov et al. Plasma Phys. Control. Fusion. (2010).
%
% USAGE: snowFinder.m
%
% INPUTS:
%
%   psizr....matrix with dimensions (nz x nr) containing the magnetic flux
%            at nz vertical by nr radial grid points
%
%   rExp.....radial coordinate   defining origin of series expansion
%
%   zExp.....vertical coordinate defining origin of series expansion
%
%   rhoExp...distance of field measurement point from expansion origin
%
%   rg.......array containing the nr radial grid points
%
%   zg.......array containing the nz vertical grid points
%
%   psibry...value of magnetic flux on the plasma boundary
%
% OUTPUTS: 
%
%   rxP.....radial   coordinate of primary null (on separatrix)
%
%   zxP.....vertical coordinate of primary null (on separatrix)
%
%   rxS.....radial   coordinate of secondary null
%
%   zxS.....vertical coordinate of secondary null
%
%   rSnow...snow centroid radial coordinate
%
%   zSnow...snow centroid vertical coordinate
%
%   dR......centroid-to-primary radial distance
%
%   dZ......centroid-to-primary vertical distance
%
%   rho.....snow radius sqrt(dR*dR + dZ*dZ)
%
%   theta...angular position of primary (0 < theta < 360deg)
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 11/16/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 11/16/2017
%
%..........................................................................

%..........................................
% Initialize variables and analyze the grid

% Grid point weighting matrix for cubic convolution (Keys, IEEE 1981)

mx = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1];

% Convert rg and zg to column vectors if necessary

if size(rg,1) == 1
    rg = rg';
end
if size(zg,1) == 1
    zg = zg';
end

dr = rg(2)-rg(1);
dz = zg(2)-zg(1);

% Convert psizr to (nz x nr) if necessary

nr = length(rg);
nz = length(zg);

if size(psizr,1) == 1 || size(psizr,2) == 1
    psizr = reshape(psizr, nz, nr);
end

%.....................................................................
% Calculate the locations of 3 sample points centered around snowflake

alpha0   = pi/6;
angleInc = (2*pi)/3;

r1 = rExp + rhoExp*cos(alpha0);
z1 = zExp + rhoExp*sin(alpha0);

r2 = rExp + rhoExp*cos(alpha0 + angleInc);
z2 = zExp + rhoExp*sin(alpha0 + angleInc);

r3 = rExp + rhoExp*cos(alpha0 + 2*angleInc);
z3 = zExp + rhoExp*sin(alpha0 + 2*angleInc);

rVec = [r1; r2; r3];
zVec = [z1; z2; z3];

scatter(rVec,zVec,'filled');

%.........................................
% Calculate Br and Bz at the sample points

bVec = zeros(6,1);

for ii = 1:3
    
    % Locate four grid points in each direction around the query point
    
    ir = find(rg < rVec(ii), 1, 'last');
    iz = find(zg < zVec(ii), 1, 'last');
    
    ii1 = [nz*(ir-2)+iz-1 nz*(ir-2)+iz nz*(ir-2)+iz+1 nz*(ir-2)+iz+2];
    ii2 = [nz*(ir-1)+iz-1 nz*(ir-1)+iz nz*(ir-1)+iz+1 nz*(ir-1)+iz+2];
    ii3 = [nz*(ir+0)+iz-1 nz*(ir+0)+iz nz*(ir+0)+iz+1 nz*(ir+0)+iz+2];
    ii4 = [nz*(ir+1)+iz-1 nz*(ir+1)+iz nz*(ir+1)+iz+1 nz*(ir+1)+iz+2]; 
    
    F = [psizr(ii1(1)) psizr(ii1(2)) psizr(ii1(3)) psizr(ii1(4)); ...
         psizr(ii2(1)) psizr(ii2(2)) psizr(ii2(3)) psizr(ii2(4)); ...
         psizr(ii3(1)) psizr(ii3(2)) psizr(ii3(3)) psizr(ii3(4)); ...
         psizr(ii4(1)) psizr(ii4(2)) psizr(ii4(3)) psizr(ii4(4))  ...
    ];

    % Normalize the r and z intervals
     
    tr = (rVec(ii) - rg(ir))/dr;
    tz = (zVec(ii) - zg(iz))/dz;
    
    % Interpolate to find dpsiDR and dpsiDZ
     
    b0 = [1 tr tr^2 tr^3]*mx*F(:,1);
    b1 = [1 tr tr^2 tr^3]*mx*F(:,2);
    b2 = [1 tr tr^2 tr^3]*mx*F(:,3);
    b3 = [1 tr tr^2 tr^3]*mx*F(:,4);
    
    b0_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,1);
    b1_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,2);
    b2_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,3);
    b3_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,4);
    
    psi_r = [1 tz tz^2 tz^3]*mx*[b0_r b1_r b2_r b3_r]';
    psi_z = ([0 1 2*tz 3*tz^2]/dz)*mx*[b0 b1 b2 b3]';
    
    bVec(ii)   = -1/(2*pi) * 1/rVec(ii) * psi_z;
    bVec(ii+3) =  1/(2*pi) * 1/rVec(ii) * psi_r;
    
end

%.............................................................
% Create arrays of normalized r and z sample point coordinates

sig1 = (r1 - rExp) / rExp;
zet1 = (z1 - zExp) / rExp;

sig2 = (r2 - rExp) / rExp;
zet2 = (z2 - zExp) / rExp;

sig3 = (r3 - rExp) / rExp;
zet3 = (z3 - zExp) / rExp;

sigVec = [sig1; sig2; sig3];
zetVec = [zet1; zet2; zet3];

%.....................................
% Solve for the expansion coefficients

m = zeros(6,6);

for ii = 1:3
    den = (1 + sigVec(ii))*power(rExp,2);
    m(ii,1)   = 0;
    m(ii,2)   = -1 / den;
    m(ii,3)   = -2 * sigVec(ii) / den;
    m(ii,4)   = -2 * zetVec(ii) / den;
    m(ii,5)   =  6 * sigVec(ii) * zetVec(ii) / den;
    m(ii,6)   =  3 * ( power(sigVec(ii),2) - power(zetVec(ii),2) ) / den;
    m(ii+3,1) =  1 / den;
    m(ii+3,2) =  0;
    m(ii+3,3) = (2 * zetVec(ii)) / den;
    m(ii+3,4) = -2 * sigVec(ii) / den;
    m(ii+3,5) =  3 * ( power(sigVec(ii),2) - power(zetVec(ii),2) ) / den;
    m(ii+3,6) = -6 * sigVec(ii) * zetVec(ii) / den;
end

x = m\bVec;

l1 = x(1); l2 = x(2); q2 = x(3); q3 = x(4); c1 = x(5); c4 = x(6);

q1 = -q3;
c2 = -3*c4;
c3 = -3*c1;

%..............................
% Solve for auxiliary variables

sigma0 = (q3*c1 + q2*c4)/(3*(c1*c1 + c4*c4));
zeta0  = (q2*c1 - q3*c4)/ (3*(c1*c1 + c4*c4));

P = (l2*c4 - l1*c1)/(3*(c1*c1 + c4*c4)) + sigma0*sigma0 - zeta0*zeta0;
Q = (l2*c1 + l1*c4)/(6*(c1*c1 + c4*c4)) + sigma0*zeta0;

sigma1p = sqrt(0.5*P + sqrt(0.25*P*P + Q*Q));
zeta1p  = sign(Q)*sqrt(-0.5*P + sqrt(0.25*P*P + Q*Q));

%..............................................................
% Compare fluxes at the two nulls to determine which is primary

sigma1 = sigma0 + sigma1p;
zeta1  = zeta0  + zeta1p;

sigma2 = sigma0 - sigma1p;
zeta2  = zeta0  - zeta1p;

psi1 = l1*sigma1 + l2*zeta1 + q1*sigma1*sigma1 + 2*q2*sigma1*zeta1    + ...
    q3*zeta1*zeta1 + c1*sigma1*sigma1*sigma1 + c2*sigma1*sigma1*zeta1 + ...
    + c3*sigma1*zeta1*zeta1 + c4*zeta1*zeta1*zeta1;

psi2 = l1*sigma2 + l2*zeta2 + q1*sigma2*sigma2 + 2*q2*sigma2*zeta2    + ...
    q3*zeta2*zeta2 + c1*sigma2*sigma2*sigma2 + c2*sigma2*sigma2*zeta2 + ...
    + c3*sigma2*zeta2*zeta2 + c4*zeta2*zeta2*zeta2;

%.................................
% Calculate the snowflake location

rx1 = rExp + rExp*sigma1;
zx1 = zExp + rExp*zeta1;

rx2 = rExp + rExp*sigma2;
zx2 = zExp + rExp*zeta2;

if nargin == 7
    
    dpsi1 = abs(psi1 - psibry);
    dpsi2 = abs(psi2 - psibry);
    
    if dpsi2 > dpsi1 % null one is primary
        
        rxP = rx1;
        zxP = zx1;
        rxS = rx2;
        zxS = zx2;
        
        dR = rExp*sigma1p;
        dZ = rExp*zeta1p;
        
    else            % null two is primary
        
        rxP = rx2;
        zxP = zx2;
        rxS = rx1;
        zxS = zx1;
        
        dR = rExp*(-sigma1p);
        dZ = rExp*(-zeta1p);
        
    end
    
else
    
    if zx1 > zx2 % null one is primary
        
        rxP = rx1;
        zxP = zx1;
        rxS = rx2;
        zxS = zx2;
        
        dR = rExp*sigma1p;
        dZ = rExp*zeta1p;
    
    else         % null two is primary
        
        rxP = rx2;
        zxP = zx2;
        rxS = rx1;
        zxS = zx1;
        
        dR = rExp*(-sigma1p);
        dZ = rExp*(-zeta1p);
        
    end
    
end

rSnow = (rxP + rxS)/2;
zSnow = (zxP + zxS)/2;

rho   = sqrt(dR*dR + dZ*dZ);
theta = atan(dR/dZ);

theta = rad2deg(theta);

end
