% Treats the inverse mapping as a cost function for optimization. 
% Designs a new eq, simulates heat flux, and evaluates cost function based
% on comparison to IR heat flux.

% opts: optional struct to use strike point (instead of x-pt) constraints 
%   opts.constrain_sp = flag to use strike pt
%   opts.sp = [rstrike1 rstrike2 zstrike1 zstrike2] pts to constrain

% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
function eq = designeq_ml_beta(xp,shot,time_ms,init,opts)

if ~exist('opts','var'), opts.constrain_sp = 0; end 

persistent iter
if isempty(iter), iter = 0; end
if iter == 0, clear init spec config gsdesign; end
iter = iter + 1;

root = '/u/jwai/d3d_snowflake_2020/current/';
[rxP,rxS,zxP,zxS] = unpack(xp);

%--------------------------
% DESIGN NEW EQUILIBRIUM
%--------------------------

% Load tokamak definition
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_6565.mat';
load([tokdir tokdata]);

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

%........................
% x-pt or strike point constraints
        
if opts.constrain_sp
  spec.targets.rx = [rxP rxS];
  spec.targets.zx = [zxP zxS];
  spec.weights.x = [1e3 1e3];  
  spec.targets.rsep(end+1:end+2) = opts.sp(1:2);
  spec.targets.zsep(end+1:end+2) = opts.sp(3:4);
  spec.weights.sep(end+1:end+2) = 1e8;
else  
  spec.locks.rx = [rxP rxS];
  spec.locks.zx = [zxP zxS];    
end

spec.targets.rbdef = rxP;
spec.targets.zbdef = zxP;
spec.weights.bdef  = 1;

%.....................
% Specify the profiles

config.constraints = 1;  % allow for scaling/peaking of profiles
config.pres0 = init.pres;
config.fpol0 = init.fpol;

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


% spec.targets.betap = 
% spec.targets.li = 
% spec.weights.betap = 1;
% spec.weights.li = 1;

%  give low weight to divertor control coils (see tok_data_struct.ccnames)
iDivertor = [11 15 16 19 20]; % F9A, F4B, F5B, F8B, F9B

spec.weights.ic(1:20) = 1;
spec.weights.ic(iDivertor) = 0.1;

%.......................
% Obtain new equilibrium
eq = gsdesign(spec, init, config);

end






















