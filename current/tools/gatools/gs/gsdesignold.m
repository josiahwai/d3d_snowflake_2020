% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   eq = gsdesign(spec, init, config)
% 
%   PURPOSE: Design a 2-D (Grad-Shafranov) equilibrium by minimizing
%            a cost function: sum(weights.param*(param-targets.param))
% 
%   INPUTS: spec, specification, a structure with:
%                 targets, weights, limits, locks
% 
%           init, initial equilibrium
% 
%         config, Toksys description of the tokamak, and options
% 
%   OUTPUTS: eq, equilibrium 
% 
% For detailed information, run gsdesign without inputs
% List of scripts that demonstrate use of gsdesign:
% gsdesign_demo_d3d_DN   - Design of EAST-like DN plasma for DIII-D
% gsdesign_demo_d3d_DSNF - Design of a double snowflake plasma
% See also gsdesign_recreate_init, gsdesign_iss
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%