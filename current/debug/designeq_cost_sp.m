% Treats the inverse mapping as a cost function for optimization. 
% Designs a new eq, simulates heat flux, and evaluates cost function based
% on comparison to IR heat flux.


% EFIT: xp = [1.1056 1.1603 -1.2072 -1.2597];
% r_spd = [1.0160    1.0664    1.0664    1.3042];
% z_spd = [-1.1078   -1.2745   -1.2745   -1.3630];


% function J = designeq_cost_sp
ccc
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

  
% ----------------------
% USER / FUNCTION INPUTS
% ----------------------

% 165288: 4200
% xp = [1.1056 1.1603 -1.2072 -1.2597];
% r_spd = [1.0160    1.0664    1.3042] + [0 0 0];
% z_spd = [-1.1078   -1.2745  -1.3630] + [0 0 0];

% 165288: 4500
% xp = [1.1221    1.1513   -1.1679   -1.2897];
% r_spd = [1.0160  1.0698  1.3037];
% z_spd = [-1.1069 -1.2780  -1.3630] + [.03 0 0];

% 155355: 4500
xp = [1.1575    1.1866   -1.1958   -1.3265];
r_spd = [1.0160    1.0685];
z_spd = [-1.0776   -1.2767];

shot = 155355;
time_ms = 4500;
[rxP,rxS,zxP,zxS] = unpack(xp);

% persistent iSim
% if isempty(iSim), iSim = 0; end
% iSim = iSim + 1;
iSim = 1;

%--------------------------
% DESIGN NEW EQUILIBRIUM
%--------------------------
fprintf('\n\nshot: %d time: %d \n\n', shot, time_ms);

% Load tokamak definition
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_6565.mat';
load([tokdir tokdata]);

% Load the equilibrium
eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
timerange = [time_ms+16 time_ms-16]/1000;
init = read_eq(shot, timerange, eqdir);

% config settings
config = tok_data_struct;
config.max_iterations = 30;
spec.max_iterations = 30;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Specify where the flux should equal the boundary flux
spec.targets.rsep = init.gdata.rbbbs([1:60 82:89]);
spec.targets.zsep = init.gdata.zbbbs([1:60 82:89]);
spec.weights.sep  = ones(length(spec.targets.rsep),1) * 1e3;

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

%........................
% x-pt prediction 

% debuggg
%================================

% spec.locks.rsep = r_spd;
% spec.locks.zsep = z_spd;
spec.targets.rsep = [spec.targets.rsep' r_spd];
spec.targets.zsep = [spec.targets.zsep' z_spd];
spec.weights.sep  = [spec.weights.sep' 1e6 1e6];

% spec.locks.rx = [rxP rxS];
% spec.locks.zx = [zxP zxS];

spec.targets.rx = [rxP rxS];
spec.targets.zx = [zxP zxS];
spec.weights.x = [1 1e5] * 1;
%================================

spec.targets.rbdef = rxP;
spec.targets.zbdef = zxP;
spec.weights.bdef  = 1;

%.....................
% Specify the profiles

config.constraints = 1;  % allow for scaling/peaking of profiles
config.pres0 = init.gdata.pres;
config.fpol0 = init.gdata.fpol;

%.................................
% Specify coil circuit connections

pp_dir = [root 'inputs/coil_data/' num2str(shot)];
pp_data = ['/PP_' num2str(shot) '.mat'];

load([pp_dir pp_data])
spec.buscode = [0 0 PP.bus_code];


%....................
% Load coil currents

coil_dir = [root 'inputs/coil_data/' num2str(shot)];
coil_data = ['/coil_currents_' num2str(shot) '.mat']; % loads currents, t
load([coil_dir coil_data]);

% smooth coil currents slightly
cc_smooth = zeros(size(coil_currents));
nc = min(size(coil_currents));

for i = 1:nc
    cc_smooth(i,:) = smooth(coil_currents(i,:), 10); 
end

[~,k] = min(abs(t-time_ms));

spec.locks.ic = nan*ones(20,1);
spec.targets.ic = cc_smooth(:,k);


%  give low weight to divertor control coils (see tok_data_struct.ccnames)
iDivertor = [11 15 16 19 20]; % F9A, F4B, F5B, F8B, F9B

spec.weights.ic(1:20) = 1;
spec.weights.ic(iDivertor) = 0.1;


%.......................
% Obtain new equilibrium
eq = gsdesign(spec, init.gdata, config);


% run heatsim
opts.iSim = 2;
opts.plotit = 1;
opts.saveit = 1;
opts.saveDir = '/u/jwai/d3d_snowflake_2020/current/debug/';
opts.root = root;

J = heatsim_cost(eq,shot,time_ms,opts);

if isnan(J) || isinf(J)
    J = 1.0;  % large penalty
end

figure(11)
plot(xp(1:2), xp(3:4), 'xr', 'Markersize', 10, 'LineWidth', 4)
r = [.1260    1.1480    1.1362    1.1401    1.1186 1.1312    1.1186    1.1312];
z = [ -1.1105   -1.2682   -1.0932   -1.2812   -1.1100 -1.2710   -1.1100   -1.2710];
plot(r, z, 'xk', 'Markersize', 10, 'LineWidth', 4)
scatter(r_spd,z_spd,'k','filled')
























