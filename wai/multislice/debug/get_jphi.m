% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
shot = 165288;
time_ms = 4000;  

plotit = 1;


% ------------------------------------------
% OBTAIN SNOWFLAKE AND HEAT FLUX PARAMETERS
% ------------------------------------------

addpath(genpath('/u/jwai/d3d_snowflake_2019_wai'));
addpath(genpath('/u/jwai/d3d_snowflake_2019'));


% Load tokamak definition
load('/u/jwai/d3d_snowflake_2019_wai/hf_constrained_eq/mat_files/d3d_obj_mks_struct_129129.mat')


% Compute total length of the limiter [m]
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
    

% Load the equilibrium
eq_dir = strcat('/u/jwai/d3d_snowflake_2019_wai/multislice/eq_unconstrained/', int2str(shot), '/cake');
eq = read_eq(shot, time_ms/1000, eq_dir);


%================================================================


jphi = reshape(eq.gdata.jphi , 129^2, 1);
psiN = (-reshape(eq.gdata.psizr, 129^2, 1) + eq.gdata.psimag) / (eq.gdata.psimag - eq.gdata.psibry);

close all
compareEdge
scatter(psiN, jphi, 'k', 'filled')










