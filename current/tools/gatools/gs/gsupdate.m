% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   To initialize: gsupdate(x, e, c) or gsupdate(x, e, c, 1)
%     x becomes initial x, e.psizr becomes initial flux
%   To update: gsupdate(x) or gsupdate(x,e,c,0)
%     uses old e, r in memory to calculate psizr for state x
%   To get all outputs: [d,e,r,b,p] = gsupdate(x)
%   To recycle any cached d for state near x: d = gsupdate(x)
% 
%   PURPOSE: Update d,e,r,b,p to state x or d to near state x
% 
%   INPUTS: x, new state vector (initialize with gsinit.m)
%           e, equilibrium 
%           c, configuration created by gsconfig.m          
% 
%   OUTPUTS: d, dynamics (to evolve x and calculate outputs)
%            e, equilibrium (common quantities)
%            r, response (how the equilibrium changes with x)
%            b, boundary (more detailed information than in e)
%            p, profiles (c.np contours from axis to boundary)
%