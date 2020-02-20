%...............................................................................
%
% USAGE: gsdesign_d3d_165286_snowMinus.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 03/13/2019
%
% PURPOSE: Design a snowflake minus equilibrium for DIII-D based on shot 165286.
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 03/13/2019
%
%...............................................................................

clear all

shot = 165286;
time = 3500;

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

efit_dirname = '/u/pvail/d3d_snowflake_2019/gsedge/';
    
init = read_eq(shot, time/1000, efit_dirname);

psimag_efit = init.gdata.psimag;
psibry_efit = init.gdata.psibry;

psi_efit = linspace(psimag_efit, psibry_efit, 129);

psibar_efit = (psi_efit - psimag_efit)./(psibry_efit - psimag_efit);

pprime_efit = init.gdata.pprime;
ffprim_efit = init.gdata.ffprim;

%......................................................
% Specify where the flux should equal the boundary flux

spec.targets.rsep = init.gdata.rbbbs([1:60 86:91]);
spec.targets.zsep = init.gdata.zbbbs([1:60 86:91]);
spec.weights.sep = ones(length(spec.targets.rsep),1);

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

spec.targets.li = 1.016;
spec.weights.li = 1;

spec.targets.betap = 1.095;
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

spec.targets.ic(3) = -2000;  % F1A
spec.targets.ic(4) = -2400;  % F2A
spec.targets.ic(5) =  0;     % F3A
spec.targets.ic(6) =  300;   % F4A
spec.targets.ic(7) =  100;   % F5A
spec.targets.ic(8) = -3450;  % F6A
spec.targets.ic(9) = -4000;  % F7A
spec.targets.ic(10) = 45;    % F8A
spec.targets.ic(11) = 0;     % F9A

spec.targets.ic(12) = -2320; % F1B
spec.targets.ic(13) = -1250; % F2B
spec.targets.ic(14) = -1670; % F3B
spec.targets.ic(15) =  3040; % F4B
spec.targets.ic(16) = -1000; % F5B
spec.targets.ic(17) = -3100; % F6B
spec.targets.ic(18) = -6000; % F7B
spec.targets.ic(19) =  2850; % F8B
spec.targets.ic(20) =  0;    % F9B

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
