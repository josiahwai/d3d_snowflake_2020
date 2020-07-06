%   USAGE:   gs_helical_voltage
% 
%   PURPOSE: Calculate voltage along B*sign(bzero) for contours
% 
%   INPUTS: etacont, parallel resistivity at points rcont, zcont
%           Create rcont, zcont and other inputs with:
%             gs_eq_analysis
%             gs_response
%             gs_trace_contours
%             gs_contour_response
%             gs_contour_profiles
%           Then create etacont before calling this script
% 
%   OUTPUTS: vresc, resistive voltage per toroidal turn
%            dvindcdxdot, inductive voltage as function of xdot
% 	    Projcurprofeqs, projection onto nkn+2 equations
% 
%   METHOD:  The contours are stationary. The outputs can be used
%            to calculate surface-integrals of dB/dt + curl E = 0
%            for surfaces bounded by field lines at fluxes psibarc,
%            using that E_para = eta_para*j_para along field lines
%