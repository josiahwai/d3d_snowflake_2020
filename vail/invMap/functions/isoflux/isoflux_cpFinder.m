function cpRZ = isoflux_cpFinder(psizr, psibry, rg, zg, segData)
%
% ISOFLUX_CPFINDER
%
%   Locate the isoflux control points for a given plasma equilibrium.
%
%   Determine the (r,z) coordinates of intersection points between the
%   plasma boundary and selected isoflux control segments.
%
% USAGE: isoflux_cpFinder.m
%
% INPUTS:
%
%   psizr....matrix with dimensions (nz x nr) containing the magnetic flux
%            at nz vertical by nr radial grid points
%
%   psibry...flux on the plasma boundary
%
%   rg.......array containing the nr radial grid points
%
%   zg.......array containing the nz vertical grid points
%
%   segData...matrix with dimensions (numSegs x 4) where each row contains
%             [rstart rend zstart zend] for a control segment
%
% OUTPUTS: 
%
%   cpRZ......matrix with dimensions (numSegs x 2) where each row contains
%             [cpR cpZ] for a control point
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 06/14/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/14/2017
%
%..........................................................................

numSegs = size(segData,1);

cpRZ = zeros(numSegs,2);

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

for ii = 1:numSegs
    
    % Compute slope of the control segment
    
    rS = segData(ii,1);
    zS = segData(ii,3);
    rE = segData(ii,2);
    zE = segData(ii,4);
    
    dR = rE - rS;
    dZ = zE - zS;
    m = dZ/dR;
    z0 = zS - m*rS;
    
    % Use center of control segment as guess for control point location
    
    cpR = rS +   dR/2; 
    cpZ = zS + m*dR/2;
    
    jj  = 20;  % maximum of 20 iterations to find the control point
    err = inf;
    
    while jj > 0 && err > 1e-10
        
        jj = jj-1;
        
        % Locate four grid points in each direction around the query point
        
        ir = find(rg < cpR, 1, 'last');
        iz = find(zg < cpZ, 1, 'last');
        
        ii1 = nz*(ir-2) + [iz-1 iz iz+1 iz+2];
        ii2 = nz*(ir-1) + [iz-1 iz iz+1 iz+2];
        ii3 = nz*(ir+0) + [iz-1 iz iz+1 iz+2];
        ii4 = nz*(ir+1) + [iz-1 iz iz+1 iz+2];
        
        F = [psizr(ii1(1)) psizr(ii1(2)) psizr(ii1(3)) psizr(ii1(4)); ...
             psizr(ii2(1)) psizr(ii2(2)) psizr(ii2(3)) psizr(ii2(4)); ...
             psizr(ii3(1)) psizr(ii3(2)) psizr(ii3(3)) psizr(ii3(4)); ...
             psizr(ii4(1)) psizr(ii4(2)) psizr(ii4(3)) psizr(ii4(4))  ...
        ];
    
        % Normalize the r and z intervals
        
        tr = (cpR - rg(ir))/dr;
        tz = (cpZ - zg(iz))/dz;
        
        % Interpolate to find psi at the control point
        
        b0 = [1 tr tr^2 tr^3]*mx*F(:,1);
        b1 = [1 tr tr^2 tr^3]*mx*F(:,2);
        b2 = [1 tr tr^2 tr^3]*mx*F(:,3);
        b3 = [1 tr tr^2 tr^3]*mx*F(:,4);
        
        b0_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,1);
        b1_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,2);
        b2_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,3);
        b3_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,4);
        
        psiCP = [1 tz tz^2 tz^3]*mx*[b0 b1 b2 b3]';
        
        psiCP_r = [1 tz tz^2 tz^3]*mx*[b0_r b1_r b2_r b3_r]';
        psiCP_z = ([0 1 2*tz 3*tz^2]/dz)*mx*[b0 b1 b2 b3]';
        
        delta = -inv([psiCP_r psiCP_z; -m 1])*[psiCP-psibry; cpZ-m*cpR-z0];
        cpR = cpR + delta(1);
        cpZ = cpZ + delta(2);
        
        cpRZ(ii,1) = cpR;
        cpRZ(ii,2) = cpZ;
        
        err = abs(psiCP - psibry);
       
    end
    
end
