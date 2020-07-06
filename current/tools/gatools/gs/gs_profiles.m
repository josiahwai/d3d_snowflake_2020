%   USAGE:   gs_profiles
% 
%   PURPOSE: Find how much of quantities are contained within flux-
%            surfaces from axis to boundary, by accumulated sums
%            of quantities within grid elements sorted w.r.t. psibarzr
% 
%   INPUTS:   profile_fit_option (default = 2)
%               0 = just fit cumsum of psibarzr-sorted grid cells
%               1 = fit equals the calculated total value at psibar = 1
%               2 = fit also equals calculated derivative at psibar = 0
%             calculate_profile_responses, a flag, default = false
% 
%             Remaining inputs are made by gs_eq_analysis:
%             psibarzr, normalized flux on the grid
%             p0, p1, p2, p3, spline function for pressure
%             Wzr, thermal energy within grid cells
%             Wth, total thermal energy
%             Vtot, total plasma volume
%             pprimezr, pprime on the grid
%             ffprimzr, ffprim on the grid
%             cpasma, total plasma current
%             Atot, total plasma area
%             fpolzr, fpol on the grid
%             torflux, total toroidal flux within plasma
% 
%   OUTPUTS: SPLINE FUNCTIONS (poly-coefficients for each spline region):
%              v0, v1, v2, v3, spline function for volume, V(psibar)
%              w0, w1, w2, w3, w4, w5, w6, thermal energy, W(psibar)
%              a0, a1, a2, a3, spline function for area, A(psibar)
%              i0, i1, i2, i3, for toroidal current, I(psibar)
%              t0, t1, t2, t3, for toroidal flux, T(psibar)
%            Values at psibar (i.e. nr values from axis to boundary):
%              V, volume contained within flux surfaces
%              W, thermal energy contained within flux surfaces
%              A, area contained within flux surfaces
%              I, toroidal current contained within flux surfaces
%              T, toroidal flux contained within flux surfaces
%              jtav, contour-averaged current density (dI/dA)
%              rhot, square-root of normalized toroidal flux, sqrt(T/T(nr))
%              qpsi, q values (dT/dpsi)
% 	
%   METHOD: The grid points inside the plasma are sorted w.r.t. psibarzr.
%           Quantities contained within flux surfaces are found by
%           accumulated sums over the sorted grid points.
%           Spline functions are fitted to these accumulated sums.
%