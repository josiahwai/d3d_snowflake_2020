% Extensive rewrite of P.J. Vail's original snowflake divertor heat flux
% analysis codes. 
% 
% Workflow for this analysis:
% 
% 1. Ensure necessary input files are available (inputs folder)
% 2. Design a new, infrared-constrained equilibrium using this script (design_eq)
%     - obtains heat flux parameters
%     - predicts xpoint locations
%     - designs new equilibrium
%     - measures and compares edge current
% 3. Use the new eq to run a heat flux simulation using the heatsim_noEich script
%
% 4. Analysis questions to answer: Do the new eqs match heat flux better?
% Are there trends in the edge current differences? How accurate is x-point
% prediction?
%
%
% Change Modification:
%       Josiah Wai, 2/26/2020, created original file
%


% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
args = importdata('args.txt');
shot = args(1);
time_ms = args(2);

plotit = 1;
saveit = 1;
simulate_heat = 1;
    
fprintf('\n\n\nshot: %d time: %d \n\n', shot, time_ms);

% ---------------------------
% OBTAIN HEAT FLUX PARAMETERS
% ---------------------------
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

% Load tokamak definition
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);


% Load the equilibrium
eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
timerange = [time_ms+16 time_ms-16]/1000;
eq = read_eq(shot, timerange, eqdir);


% Load heat flux data q(s,t), s=distance along limiter, and t=time 
qperp_dir  = [root 'inputs/qperp/' num2str(shot) '/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];

load([qperp_dir qperp_data])  % loads q, s, and t

[~,k] = min(abs(t-time_ms));
qperp = qperp(k,:)';


% Remove the gap from s (distance along limiter)
idxGap = find(s < 170);
gap1 = s(idxGap(end));
gap2 = s(idxGap(end)+1);
dgap = gap2 - gap1;
s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;


% Index data for each SP
idx_SP1 = find(s < 115);
idx_SP3 = find(s > 145);

idx_SP2 = setdiff(1:length(s), [idx_SP1 idx_SP3]);

qperp_SP1 = qperp(idx_SP1);
qperp_SP2 = qperp(idx_SP2);
qperp_SP3 = qperp(idx_SP3);

s_SP1 = s(idx_SP1)';
s_SP2 = s(idx_SP2)';
s_SP3 = s(idx_SP3)';


% Determine magnitude and location of heat flux peak at each SP
[qmax_SP1_temp, idx_SP1] = max(qperp_SP1);
[qmax_SP2_temp, idx_SP2] = max(qperp_SP2);
[qmax_SP3_temp, idx_SP3] = max(qperp_SP3);

sDiv_qmax_SP1_temp = s_SP1(idx_SP1);
sDiv_qmax_SP2_temp = s_SP2(idx_SP2);
sDiv_qmax_SP3_temp = s_SP3(idx_SP3);


% Determine FWHM at each SP
idx11 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'first');
idx21 = find(qperp_SP1 >= qmax_SP1_temp/2, 1, 'last');

FWHM_SP1_temp = abs(s_SP1(idx11) - s_SP1(idx21));

idx12 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'first');
idx22 = find(qperp_SP2 >= qmax_SP2_temp/2, 1, 'last');

FWHM_SP2_temp = abs(s_SP2(idx12) - s_SP2(idx22));

idx13 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'first');
idx23 = find(qperp_SP3 >= qmax_SP3_temp/2, 1, 'last');

FWHM_SP3_temp = abs(s_SP3(idx13) - s_SP3(idx23));

sDiv_qmax_SP1_temp = 0.01*sDiv_qmax_SP1_temp;
sDiv_qmax_SP2_temp = 0.01*sDiv_qmax_SP2_temp;
sDiv_qmax_SP3_temp = 0.01*sDiv_qmax_SP3_temp;

FWHM_SP1_temp = 0.01*FWHM_SP1_temp;
FWHM_SP2_temp = 0.01*FWHM_SP2_temp;
FWHM_SP3_temp = 0.01*FWHM_SP3_temp;

qmax_SP1_temp = qmax_SP1_temp*(1e-6)*(100)^2;
qmax_SP2_temp = qmax_SP2_temp*(1e-6)*(100)^2;
qmax_SP3_temp = qmax_SP3_temp*(1e-6)*(100)^2;


%--------------------------
% PREDICT X-POINT LOCATIONS
%--------------------------

% train the regression tree
% HeatFluxMapping_RegressionTree_v6    

% load the regression tree
rtree_dir = [root 'design_eq/invMap/tree/v6/'];
rtree = {'rtree_rxP.mat', 'rtree_rxS.mat', 'rtree_zxP.mat', 'rtree_zxS.mat'};

for i = 1:length(rtree)
    load([rtree_dir rtree{i}])
end


% predict x-point positions
rxP_PRED = predict(rtree_rxP, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp]);
rxS_PRED = predict(rtree_rxS, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp]);
zxP_PRED = predict(rtree_zxP, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp]);
zxS_PRED = predict(rtree_zxS, [qmax_SP1_temp qmax_SP2_temp qmax_SP3_temp   ...
    sDiv_qmax_SP1_temp sDiv_qmax_SP2_temp sDiv_qmax_SP3_temp FWHM_SP1_temp ...
    FWHM_SP2_temp FWHM_SP3_temp]);

if plotit
    figure(11)
    plot(rxP_PRED, zxP_PRED, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
        'LineWidth', 3)
    hold on
    plot(rxS_PRED, zxS_PRED, 'x', 'color', [0.91 0.41 0.17], 'Markersize', 12, ...
        'LineWidth', 3)
end


%--------------------------
% DESIGN NEW EQUILIBRIUM
%--------------------------
clear spec config

objs_dir  = [root 'inputs/tok_data/'];
objs_name = 'd3d_obj_mks_struct_6565.mat';
load([objs_dir objs_name])

config = tok_data_struct;
config.max_iterations = 2;
spec.max_iterations = 5;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;


% save EFIT data
init = eq;

psizr_efit  = init.gdata.psibry;
psimag_efit = init.gdata.psimag;
psibry_efit = init.gdata.psibry;
psi_efit    = linspace(psimag_efit, psibry_efit, 129);
psibar_efit = (psi_efit - psimag_efit)./(psibry_efit - psimag_efit);
pprime_efit = init.gdata.pprime;
ffprim_efit = init.gdata.ffprim;


% Specify where the flux should equal the boundary flux
spec.targets.rsep = init.gdata.rbbbs([1:60 82:89]);
spec.targets.zsep = init.gdata.zbbbs([1:60 82:89]);
spec.weights.sep  = ones(length(spec.targets.rsep),1);


%..........................
% Specify scalar quantities
spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

% spec.targets.li = 1.015;
% spec.weights.li = 1;
% 
% spec.targets.betap = 0.932;
% spec.weights.betap = 1;


%...........................
% Specify the null locations

spec.targets.rx = [rxP_PRED rxS_PRED];
spec.targets.zx = [zxP_PRED zxS_PRED];

spec.weights.x = [1e6 1e6];

spec.targets.rbdef = rxP_PRED;
spec.targets.zbdef = zxP_PRED;
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


%........................
% Load coil currents

coil_dir = [root 'inputs/coil_data/' num2str(shot)];
coil_data = ['/coil_currents_' num2str(shot) '.mat']; % loads currents, t

load([coil_dir coil_data]);

% smooth coil currents (very noisy w.r.t time)
cc_smooth = zeros(size(coil_currents));
nc = min(size(coil_currents));

for i = 1:nc
    cc_smooth(i,:) = smooth(coil_currents(i,:), 80);
end

[~,k] = min(abs(t-time_ms));
cc0 = cc_smooth(:,k);

spec.locks.ic = nan*ones(20,1);

iA = 1:2:17; % index for A coils e.g. F1A, F2A, ..., F9A
iB = 2:2:18; % index for B coils
spec.targets.ic(3:11)  = cc0(iA);
spec.targets.ic(12:20) = cc0(iB);

%  give low weight to divertor control coils (see tok_data_struct.ccnames)
iDivertor = [11 15 16 19 20]; % F9A, F4B, F5B, F8B, F9B

spec.weights.ic(1:20) = 1;
spec.weights.ic(iDivertor) = 1e-6;


%.............
% Obtain new equilibrium

eq = gsdesign(spec, init.gdata, config);



if saveit
    out_dir = [root 'outputs/eqs/' num2str(shot) '/'];
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    out_name = ['eq_c_' num2str(shot) '_' num2str(time_ms)];
    save([out_dir out_name], 'eq');
end


% --------------------
% COMPARE EDGE CURRENT
% --------------------

jpar = eq.jpar;
psin = eq.psibar;

if plotit    
    % load previous edge current
    jdir = [root 'inputs/edge_current/' num2str(shot)];
    jdata = ['/jdata_' num2str(shot) '.mat'];
    load([jdir jdata])
    
    [~,k] = min(abs(times-time_ms));
    jpar0 = jtot(k,:);
    psin0 = psi_n;
    
    figure(33)
    hold on
    plot(psin,  jpar,  'b', 'linewidth', 2)
    plot(psin0, jpar0, 'r', 'linewidth', 2)
end


% --------------------
% SIMULATE HEAT FLUX
% --------------------
if simulate_heat
    opts.plotit = plotit;
    opts.saveit = saveit;
    opts.root = root;
    
    try
        % heat sim for the unconstrained eq (from cake)
        opts.saveDir = [root 'outputs/hfsims/cake_unconstrained/' num2str(shot) '/'];
        if ~exist(opts.saveDir, 'dir')
            mkdir(opts.saveDir);
        end
        heatsim_noEich(init.gdata, shot, time_ms, opts);        
    catch
        fprintf('Warning: CAKE unconstrained heat sim did not run. %d %dms\n', shot, time_ms);        
    end
    
    
    try
        % save the new constrained eq (based off cake)
        opts.saveDir = [root 'outputs/hfsims/cake_constrained/' num2str(shot) '/'];
        if ~exist(opts.saveDir, 'dir')
            mkdir(opts.saveDir);
        end
        heatsim_noEich(eq, shot, time_ms, opts);
    catch
        fprintf('Warning: CAKE constrained heat sim did not run. %d %dms\n', shot, time_ms);        
    end
     
    
    try
        % heat sim from the unconstrained efit01 eq
        efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
        efit_eq = read_eq(shot, time_ms/1000, efit_dir);
        opts.saveDir = [root 'outputs/hfsims/efit_unconstrained/' num2str(shot) '/'];
        if ~exist(opts.saveDir, 'dir')
            mkdir(opts.saveDir);
        end
        heatsim_noEich(efit_eq.gdata, shot, time_ms, opts);
    catch
        fprintf('Warning: EFIT unconstrained heat sim did not run. %d %dms\n', shot, time_ms);        
    end
    
end


























