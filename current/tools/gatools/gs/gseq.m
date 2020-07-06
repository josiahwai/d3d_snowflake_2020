% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [y, eq] = gseq(xc, init, config)
% 
%   PURPOSE: Solve the 2-D (free-boundary) Grad-Shafranov equation
%            CONFIGURE & INITIALIZE: xc0 = gseq([], init, config)
%            SOLVE for equilibrium: [y, eq] = gseq(xc)
% 
%   INPUTS:     xc, equilibrium states
%             init, initial equilibrium, for further instructions type:
%                   edit gseq_init.m
%           config, tokamak description and options, to set these type:
%                   edit gseq_config.m
% 
%   OUTPUTS: y, outputs defined in config.outputs, 
%               gseq('index_in_y') returns indices for named outputs
%           eq, equilibrium in TokSys format
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%