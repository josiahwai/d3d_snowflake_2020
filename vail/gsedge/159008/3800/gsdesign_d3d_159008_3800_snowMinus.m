%...............................................................................
%
% USAGE: gsdesign_d3d_159008_3800_snowMinus.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 03/25/2019
%
% PURPOSE: Design a snowflake minus equilibrium for DIII-D based on shot 159008.
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 03/25/2019
%
%...............................................................................

clear all

shot = 159008;
time = 3800;

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

efit_dirname = '/u/pvail/d3d_snowflake_2019/gsedge/159008/3800/';
    
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

% spec.targets.li = 1.051;
% spec.weights.li = 1;
% 
% spec.targets.betap = 0.907;
% spec.weights.betap = 1;

% Specify the profiles

config.pres0 = init.gdata.pres;
 
config.fpol0 = init.gdata.fpol;

%.................................
% Specify coil circuit connections

pp_dir  = '/u/pvail/d3d_snowflake_2019/gsedge/159008/';
pp_name = 'ppconfig_159008.mat';

load([pp_dir pp_name])

spec.buscode = [0 0 ppconfig.bus_code];

%........................
% Define PF coil currents

spec.locks.ic = nan*ones(29,1);

spec.targets.ic(3)  = -2007;  % F1A
spec.targets.ic(4)  = -2172;  % F2A
spec.targets.ic(5)  =  44;    % F3A
spec.targets.ic(6)  =  276;   % F4A
spec.targets.ic(7)  =  97;    % F5A
spec.targets.ic(8)  = -3418;  % F6A
spec.targets.ic(9)  = -3559;  % F7A
spec.targets.ic(10) =  161;   % F8A
spec.targets.ic(11) = -51;    % F9A

spec.targets.ic(12) = -2053;  % F1B
spec.targets.ic(13) = -1568;  % F2B
spec.targets.ic(14) = -2020;  % F3B
spec.targets.ic(15) =  3514;  % F4B
spec.targets.ic(16) = -1105;  % F5B
spec.targets.ic(17) = -2966;  % F6B
spec.targets.ic(18) = -5744;  % F7B
spec.targets.ic(19) =  2994;  % F8B
spec.targets.ic(20) =  50;    % F9B

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

legend('EFIT', 'gsdesign', 'Location', 'best')

%.............
subplot(2,1,2)
plot(psibar_efit, ffprim_efit, '-r', 'LineWidth', 2)
hold on
plot(eq.psibar, eq.ffprim, '-b', 'LineWidth', 2)

xlabel('Normalized Poloidal Flux \psi_N')
title('ffprim')
grid on
