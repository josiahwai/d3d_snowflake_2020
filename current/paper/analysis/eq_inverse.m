% Change in x-pts causes change in jboot. Test the inverse relationship. If
% jboot is given, is the new x-pt position as predicted?
ccc

shot = 155354;
time_ms = 3727;

load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/3727/eqs.mat')
psin = eqs{2}.psibar;

figure
subplot(3,1,1)
hold on
plot(psin, eqs{2}.jpar, 'b')
plot(psin, eqs{3}.jpar, 'r')
subplot(3,1,2)
hold on
plot(psin, eqs{2}.pres, 'b')
plot(psin, eqs{3}.pres, 'r')
subplot(3,1,3)
hold on
plot(psin, eqs{2}.fpol, 'b')
plot(psin, eqs{3}.fpol, 'r')

% profiles
fpol_t0 = eqs{2}.fpol;
pres_t0 = eqs{2}.pres;
fpol_tf = eqs{3}.fpol;
pres_tf = eqs{3}.pres;

root = '/u/jwai/d3d_snowflake_2020/current/';

% ==========================
% DESIGN NEW EQUILIBRIUM
% ==========================
clear spec init config

% Load tokamak definition
tokdata = 'd3d_obj_mks_struct_6565.mat';
load(tokdata);

% Load an initial equilibrium, if not supplied
eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
timerange = [time_ms+5 time_ms-5]/1000;
if ~exist('init','var')
  init = read_eq(shot, timerange, eqdir);
  init = init.gdata;
end

% config settings
config = tok_data_struct;
config.max_iterations = 30;
spec.max_iterations = 30;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Specify where the flux should equal the boundary flux
rbbbs = init.rbbbs;
zbbbs = init.zbbbs;
k = find(zbbbs > -0.65 | rbbbs > 1.85);
spec.targets.rsep = rbbbs(k);
spec.targets.zsep = zbbbs(k);
spec.weights.sep  = ones(length(spec.targets.rsep),1) * 1e3;

%..........................
% Specify scalar quantities

spec.targets.cpasma = init.cpasma;
spec.weights.cpasma = 1;
    
%.....................
% Specify the profiles

config.constraints = 1; 
config.pres0 = pres_tf';
config.fpol0 = fpol_tf';

% config.nkn = -2;

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
eq = gsdesign(spec, init, config);



















