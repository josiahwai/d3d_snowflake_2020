% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [spec,config] = gsdesign_recreate_init(init,tok_data_struct)
%            [spec,config,eq] = gsdesign_recreate_init(init,tok_data_struct)
% 
%   PURPOSE: Create inputs spec and config that make
%            init == gsdesign(spec, init, config)
% 
%   INPUTS:  init, an equilibrium to recreate
%            tok_data_struct, TokSys description of the tokamak
% 
%   OUTPUTS: spec, a specification of targets to use as input to gsdesign
%            config, configuration data to use as third input to gsdesign
%            eq = gsdesign(spec,init,config), equilibrium resembling init
% 
%   RESTRICTIONS: spec does NOT specify circuit constraints
%                 run gsdesign without inputs for full documentation
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%