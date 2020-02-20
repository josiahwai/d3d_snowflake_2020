%...............................................................................
%
% USAGE: gsdesign_d3d_165288_4000_EQHF.m
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

clear spec config

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

topdir = '/u/pvail/d3d_snowflake_2019/heatfluxconstrain/165288/4000/';

efit_dirname = [topdir 'EFIT01'];
    
init = read_eq(shot, time/1000, efit_dirname);

psizr_efit = init.gdata.psibry;

psimag_efit = init.gdata.psimag;
psibry_efit = init.gdata.psibry;

psi_efit = linspace(psimag_efit, psibry_efit, 129);

psibar_efit = (psi_efit - psimag_efit)./(psibry_efit - psimag_efit);

pprime_efit = init.gdata.pprime;
ffprim_efit = init.gdata.ffprim;

%......................................................
% Specify where the flux should equal the boundary flux

spec.targets.rsep = init.gdata.rbbbs([1:60 86:90]);
spec.targets.zsep = init.gdata.zbbbs([1:60 86:90]);
spec.weights.sep = 1e3*ones(length(spec.targets.rsep),1);

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

spec.targets.li = 1.015;
spec.weights.li = 1;

spec.targets.betap = 0.932;
spec.weights.betap = 1;

%...........................
% Specify the null locations

% rxP_PRED =  1.1302 + 0.00;
% zxP_PRED = -1.1142 + 0.05;
% 
% rxS_PRED =  1.1838 - 0.04;
% zxS_PRED = -1.2996 + 0.02;

% rxP_PRED =  1.1302 - 0.01;
% zxP_PRED = -1.1142 + 0.05;
% 
% rxS_PRED =  1.1838 - 0.04;
% zxS_PRED = -1.2996 + 0.02;

spec.targets.rx = [rxP_PRED rxS_PRED];
spec.targets.zx = [zxP_PRED zxS_PRED];

spec.weights.x = [1e6 1e6];

spec.targets.rbdef = rxP_PRED;
spec.targets.zbdef = zxP_PRED;
spec.targets.bdef  = 1;

%.....................
% Specify the profiles

config.pres0 = init.gdata.pres;
 
config.fpol0 = init.gdata.fpol;

%.................................
% Specify coil circuit connections

pp_dir  = '/u/pvail/d3d_snowflake_2019/heatfluxconstrain/165288/4000/';
pp_name = 'ppconfig_165286.mat';

load([pp_dir pp_name])

spec.buscode = [0 0 ppconfig.bus_code];

%........................
% Define PF coil currents

spec.locks.ic = nan*ones(20,1);

spec.targets.ic(3) =  -2079;  % F1A
spec.targets.ic(4) =  -2589;  % F2A
spec.targets.ic(5) =  -20;    % F3A
spec.targets.ic(6) =   283;   % F4A
spec.targets.ic(7) =   117;   % F5A
spec.targets.ic(8) =  -3749;  % F6A
spec.targets.ic(9) =  -3735;  % F7A
spec.targets.ic(10) =  37;    % F8A
spec.targets.ic(11) = -329;   % F9A

spec.targets.ic(12) = -2419; % F1B
spec.targets.ic(13) = -1192; % F2B
spec.targets.ic(14) = -1803; % F3B
spec.targets.ic(15) =  4083; % F4B
spec.targets.ic(16) = -2878; % F5B
spec.targets.ic(17) = -2092; % F6B
spec.targets.ic(18) = -7289; % F7B
spec.targets.ic(19) =  4286; % F8B
spec.targets.ic(20) = -22;   % F9B

spec.weights.ic = [1 1 1 1 1 1 1 1 1 1 1e-6 1 1 1 1e-6 1e-6 1 1 1e-6 1e-6];

%.............
% Run gsdesign

eq = gsdesign(spec, init.gdata, config);

% pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,129));
% ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,129));
% jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,129));
% 
% psibar_gs = linspace(0,1,129);

%..................
% Plot the profiles

% figure()
% 
% subplot(3,1,1)
% plot(psibar_efit, pprime_efit, '-b', 'LineWidth', 2)
% hold on
% plot(psibar_gs, pprime_gs, '-r', 'LineWidth', 2)
% 
% grid on
% xlabel('Normalized Flux \psi_N')
% title('pprime')
%     
% subplot(3,1,2)
% plot(psibar_efit, ffprim_efit, '-b', 'LineWidth', 2)
% hold on
% plot(psibar_gs, ffprim_gs, '-r', 'LineWidth', 2)
% 
% grid on
% xlabel('Normalized Flux \psi_N')
% title('ffprim')
%     
% subplot(3,1,3)
% plot(psibar_gs, jpar_gs/1e6, '-r', 'LineWidth', 2)
% 
% grid on
% xlabel('Normalized Flux \psi_N')
% title('parallel current density [MA/m^2]')
    
% save('eq_165288_4000_EQHF.mat', 'eq')
