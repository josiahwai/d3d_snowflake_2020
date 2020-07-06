% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [xsdot, y] = gsevolve(xs, u, init, config)
% 
%   PURPOSE: Simulate 2-D (Grad-Shafranov) equilibrium evolution
%            CONFIGURE & INITIALIZE: xs0 = gsevolve([], [], init, config)
%            EVOLVE the equilibrium: [xsdot, y] = gsevolve(xs, u)
% 
%   INPUTS:     xs, equilibrium states
%                u, inputs, format depends on config
%             init, initial equilibrium, for further instructions type:
%                   edit gsevolve_init.m
%           config, tokamak description and options, to set these type:
%                   edit gsevolve_config.m
% 
%   OUTPUTS: xsdot, time derivative of states (dxs/dt) if xs supplied in call
%                   OR initial states xs0 if xs not supplied in call
%                y, outputs defined in config.outputs, 
%                   gsevolve('index_in_y') returns indices for named outputs
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%