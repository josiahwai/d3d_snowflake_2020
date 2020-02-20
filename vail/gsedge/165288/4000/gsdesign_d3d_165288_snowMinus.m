%...............................................................................
%
% USAGE: gsdesign_d3d_165288_snowMinus.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 03/13/2019
%
% PURPOSE: Design a snowflake minus equilibrium for DIII-D based on shot 165288.
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 03/13/2019
%
%...............................................................................

clear all

shot = 165288;
time = 4000;

%.............................................
% Load DIII-D vacuum objects (tok_data_struct)

objs_dir  = '/u/pvail/d3d_snowflake_2019/invMap/';
objs_name = 'd3d_obj_mks_struct_6565.mat';

load([objs_dir objs_name])

config = tok_data_struct;
config.max_iterations = 30;

config.constraints = 1;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

%...........................
% Load EFIT data from gEQDSK

efit_dirname = '/u/pvail/d3d_snowflake_2019/gsedge/165288/4000/';
    
init = read_eq(shot, time/1000, efit_dirname);

psimag_efit = init.gdata.psimag;
psibry_efit = init.gdata.psibry;

psi_efit = linspace(psimag_efit, psibry_efit, 129);

psibar_efit = (psi_efit - psimag_efit)./(psibry_efit - psimag_efit);

pprime_efit = init.gdata.pprime;
ffprim_efit = init.gdata.ffprim;

%......................................................
% Specify where the flux should equal the boundary flux

spec.targets.rsep = init.gdata.rbbbs([1:60 86:89]);
spec.targets.zsep = init.gdata.zbbbs([1:60 86:89]);
spec.weights.sep = ones(length(spec.targets.rsep),1);

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

spec.targets.li = 1.015;
spec.weights.li = 1;

spec.targets.betap = 0.943;
spec.weights.betap = 1;

% Specify the profiles

config.pres0 = init.gdata.pres;
 
config.fpol0 = init.gdata.fpol;

%.................................
% Specify coil circuit connections

pp_dir  = '/u/pvail/d3d_snowflake_2019/gsedge/';
pp_name = 'ppconfig_165286.mat';

load([pp_dir pp_name])

spec.buscode = [0 0 ppconfig.bus_code];

%........................
% Define PF coil currents

spec.locks.ic = nan*ones(29,1);

% spec.targets.ic(3)  = -2040; % F1A
% spec.targets.ic(4)  = -2496; % F2A
% spec.targets.ic(5)  = -25;   % F3A
% spec.targets.ic(6)  =  281;  % F4A
% spec.targets.ic(7)  =  140;  % F5A
% spec.targets.ic(8)  = -3759; % F6A
% spec.targets.ic(9)  = -3618; % F7A
% spec.targets.ic(10) =  26;   % F8A
% spec.targets.ic(11) = -340;  % F9A
% 
% spec.targets.ic(12) = -2409; % F1B
% spec.targets.ic(13) = -1124; % F2B
% spec.targets.ic(14) = -1809; % F3B
% spec.targets.ic(15) =  3967; % F4B
% spec.targets.ic(16) = -2749; % F5B
% spec.targets.ic(17) = -2124; % F6B
% spec.targets.ic(18) = -7168; % F7B
% spec.targets.ic(19) =  1189; % F8B
% spec.targets.ic(20) =  10;   % F9B

%.............
% Run gsdesign

eq = gsdesign(spec, init.gdata, config);

%.....................
% Compare the profiles

figure(2)

%.............
subplot(2,1,1)
plot(psibar_efit, pprime_efit, '-r', 'LineWidth', 2)
hold on
plot(eq.psibar, eq.pprime, '-b', 'LineWidth', 2)

xlabel('Normalized Poloidal Flux \psi_N')
title('pprime')
grid on

legend('EFIT', 'gsdesign')

%.............
subplot(2,1,2)
plot(psibar_efit, ffprim_efit, '-r', 'LineWidth', 2)
hold on
plot(eq.psibar, eq.ffprim, '-b', 'LineWidth', 2)

xlabel('Normalized Poloidal Flux \psi_N')
title('ffprim')
grid on
