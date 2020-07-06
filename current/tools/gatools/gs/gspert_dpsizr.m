%   USAGE: [dpsizr, dpsimag, dpsibry] = gspert_dpsizr(response,perturbation,tok_data_struct,eq0,idoplot,options)
% 
%   PURPOSE: Return the perturbed flux on the grid, dpsizr,
%            given perturbation of set of independent variables defined in gspert
% 
%   INPUTS:
%     response = the plasma response calculated by gspert
%     perturbation = row vector of perturbations of the independent variables in gspert
%       May contain several columns with different perturbations 
%       If the gspert response was calculated with iconstraints = 1
%         the columns should contain [dis; dbetap; dli; dip]
%       If the gspert response was calculated with iconstraints = 2
%         the columns should contain [dis; dw; dli; dip]
%       If the gspert response was calculated with iconstraints = 3
%         the columns should contain [dis; dW; dI]
%       Here, 'is' is conductor currents (coils and vessel)
%         Should be compatible with tok_data_struct so that dpsizr_app = [mpc mpv]*dis
%       w is total thermal energy
%       W is thermal energy in a few profile points (help gspert for more info)
%       I is toroidal current in a few profile points (help gspert for more info)
%     tok_data_struct is standard toksys description of tokamak, including mutual inductances
%       should be the same that was used in gspert to obtain response
%     eq0 is the (unperturbed) equilibrium that response was calculated from
%     idoplots = (optional) if scalar, show contours of eq.psizr+dpsizr for idoplot seconds: 
%        default is 0 = don't plot
%        idoplots can be an array to plot only select time samples
%     options = optional options:
%       options.trace: if ~0, boundary is traced to make small correction of predicted dpsibry.
%         This feature is only used in the contour plot
% 
%   OUTPUTS:
%     dpsizr = the perturbed flux
%       If perturbation is only one time sample then dpsizr is a (nz,nr) matrix
%       If perturbation is matrix with nt time samples, dpsizr is a (nr*nz,nt) matrix
% 
%   RESTRICTIONS:
% 
%   METHOD:
%