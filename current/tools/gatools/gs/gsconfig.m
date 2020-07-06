% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   c = gsconfig(tok,opt)
% 
%   PURPOSE: Configure parameters used by gseq, gsdesign, gsupdate, GSevolve
% 
%   INPUTS: tok, structure with TokSys tokamak info, tok_data_struct
%           opt, options
%           For complete information type   edit gsconfig_template.m
% 
%   OUTPUTS: c, configuration data for gs codes
%            (c.info = description of variables in c)
%