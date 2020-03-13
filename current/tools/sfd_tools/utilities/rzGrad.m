function  [dFdr, dFdz] = rzGrad(F, rg, zg)
%
% RZGRAD
%
%   Script to compute the gradient of a function defined on the (r,z)
%   plasma grid.
%
% USAGE: rzGrad.m
%
% INPUTS:
%
%   F.........matrix with dimensions (nz x nr) containing the function vals
%             at nz vertical by nr radial grid points
%
%   rg........array containing the nr radial grid points
%
%   zg........array containing the nz vertical grid points
%
% OUTPUTS: 
%
%   dFdr......matrix with dimensions (nz x nr) containing the derivatives
%             of the function F at each grid point with respect to the
%             radial coordinate
%
%   dFdz......matrix with dimensions (nz x nr) containing the derivatives
%             of the function F at each grid point with respect to the
%             vertical coordinate
%
% METHOD: Computes the derivatives using centered differences for the
%         interior points and one-sided differences for the end points.
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 06/13/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/13/2017
%
%...............................................................................


%............................
% Analyze the grid parameters

nr = length(rg);
nz = length(zg);

dr = rg(2) - rg(1);
dz = zg(2) - zg(1);

%....................................
% Convert F to (nz x nr) if necessary

if size(F,1) == 1 || size(F,2) == 1
    F = reshape(F, nz, nr);
end

%.....................
% Compute the gradient

dFdr = zeros(nz,nr);
dFdz = zeros(nz,nr);

for ii = 1:nz
    for jj = 1:nr
        
        % Forward difference
        if jj == 1
            dFdr(ii,jj) = (F(ii,jj+1) - F(ii,jj))/dr; 
        % Backward difference
        elseif jj == nr
            dFdr(ii,jj) = (F(ii,jj) - F(ii,jj-1))/dr;
        % Centered difference
        else
            dFdr(ii,jj) = (F(ii,jj+1) - F(ii,jj-1))/(2*dr);
        end
        
    end
end

for jj = 1:nr
    for ii = 1:nz
        
        % Forward difference
        if ii == 1
            dFdz(ii,jj) = (F(ii+1,jj) - F(ii,jj))/dz; 
        % Backward difference
        elseif ii == nz
            dFdz(ii,jj) = (F(ii,jj) - F(ii-1,jj))/dz;
        % Centered difference
        else
            dFdz(ii,jj) = (F(ii+1,jj) - F(ii-1,jj))/(2*dz);
        end
        
    end
end
     
end
