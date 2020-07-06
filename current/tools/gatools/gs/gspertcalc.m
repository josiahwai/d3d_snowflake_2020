%   USAGE: [dpsizr,djphi,dbr,dbz,dpsizrerror] = gspertcalc(djphidpsizr,dpsizr_app,rg,zg,mpp,method,options)
% 
%   PURPOSE: Calculate plasma response to applied flux for a given djphi/dpsizr
% 
%   INPUTS:  djphidpsizr, how jphi responds to local flux perturbation [MA/m2/Wb]
%             dpsizr_app, perturbation(s) of applied flux [Wb] (size [nz,nr] or [nz*nr,:], e.g. mpc)
%                 rg, zg, vectors of R, Z coordinates for the grid [m]
%                    mpp, (optional) supply (nz*nr,nr) mutuals [H] if available.
%                 method, (optional, otherwise set automatically) calculation method (1, 2 or 3)
%                         1 = Exact solution optimized for small grids and few dpsizr_app
%                         2 = Exact solution optimized for larger grids and several dpsizr_app
%                         3 = Find solution for coarser grid and remove the error that comes from
%                             interpolating back to original grid by iterations that remove the
%                             error and associated plasma response. This method may fail if the
%                             coarse grid can't well represent the gradients in djphidpsizr
%                 options, (with method 3), structure with optional fields:
%                          xr, rg for coarse grid is rg(1:xr:nr), defaults for xr, xz depend on gradients
%                          xz, zg for coarse grid is zg(1:xz:nz), in djphidpsizr and available memory
%                          maxerror (default = 1e-12), maximum acceptable value of dpsizrerror
%                          maxiter (default = 9), maximum iterations (even if dpsizrerror > maxerror)
%                          idoplot (default = 0), plot the remaining dpsizr errors after each iteration
% 
%   OUTPUTS: dpsizr, perturbed flux [Wb] (= gscalc(djphi,rg,zg)+dpsizr_app)
%             djphi, perturbed current density [MA] (= djphidpsizr.*dpsizr)
%               dbr, perturbed radial field [T]
%               dbz, perturbed vertical field [T]
%       dpsizrerror, max(max(abs(gscalc(djphi,rg,zg)+dpsizr_app-dpsizr))), only recommended for method 3 
% 
%   METHOD: Finite differences are solved for interior and mutuals are used for edge of grid
%