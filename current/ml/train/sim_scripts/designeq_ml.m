% Treats the inverse mapping as a cost function for optimization. 
% Designs a new eq, simulates heat flux, and evaluates cost function based
% on comparison to IR heat flux.

% use opts to specify whether strike pts should be locked, not x-pts using:
% opts.lock_strike = 1; 
% opts.rstrike = [rstrike1 rstrike2]; 
% opts.zstrike = [zstrike1 zstrike2];


% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
function eq = designeq_ml(xp,shot,time_ms)
clear spec init config gsdesign

root = '/u/jwai/d3d_snowflake_2020/current/';
[rxP,rxS,zxP,zxS] = unpack(xp);

%--------------------------
% DESIGN NEW EQUILIBRIUM
%--------------------------

% Load tokamak definition
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_6565.mat';
load([tokdir tokdata]);

% Load the equilibrium
eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
timerange = [time_ms+5 time_ms-5]/1000;
init = read_eq(shot, timerange, eqdir);

% config settings
config = tok_data_struct;
config.max_iterations = 30;
spec.max_iterations = 30;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Specify where the flux should equal the boundary flux
rbbbs = init.gdata.rbbbs;
zbbbs = init.gdata.zbbbs;
k = find(zbbbs > -0.65 | rbbbs > 1.85);
spec.targets.rsep = rbbbs(k);
spec.targets.zsep = zbbbs(k);
spec.weights.sep  = ones(length(spec.targets.rsep),1) * 1e3;

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

%........................
% x-pt prediction
spec.locks.rx = [rxP rxS];
spec.locks.zx = [zxP zxS];

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
























