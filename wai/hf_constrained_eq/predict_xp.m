% Load the heat flux data, convert to regression tree predictors, and
% predict xpoint locations

shot = 165288;
time_ms = 4000;

% load qperp data for the shot (qperp, t, s)
load(['/u/jwai/d3d_snowflake_2019_wai/qperp/qperp_' int2str(shot) '.mat']);

[~,k] = min(abs(t-time_ms));
qperp = qperp_all(k,:)';


% Remove the gap
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


% Determine Regression Tree
HeatFluxMapping_RegressionTree_v6


% Scale the data to match the total power used in heat flux simulations
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


%........................................
% Design the constrained CAKE equilibrium
%........................................
clear spec config

% Load DIII-D vacuum objects (tok_data_struct)
tok_dir = '/u/jwai/d3d_snowflake_2019_wai/hf_constrained_eq/mat_files/';
obj = 'd3d_obj_mks_struct_6565.mat';
load([tok_dir obj])

config = tok_data_struct;
config.max_iterations = 30;
config.constraints = 1;

config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;


% Load the (unconstrained) CAKE equilibrium
cake_dir = '/u/jwai/d3d_snowflake_2019_wai/gfiles/165288/from-anthony';
%cake_dir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/EFIT01';
%cake_dir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/CAKEHF';

init = read_eq(shot, time_ms/1000, cake_dir);

psizr_efit = init.gdata.psibry;

psimag_efit = init.gdata.psimag;
psibry_efit = init.gdata.psibry;

psi_efit = linspace(psimag_efit, psibry_efit, 129);

psibar_efit = (psi_efit - psimag_efit)./(psibry_efit - psimag_efit);




pprime_efit = init.gdata.pprime;
ffprim_efit = init.gdata.ffprim;


% Specify boundary targets

spec.targets.rsep = init.gdata.rbbbs([1:60 82:89]);
spec.targets.zsep = init.gdata.zbbbs([1:60 82:89]);
spec.weights.sep = ones(length(spec.targets.rsep),1);


% Specify scalar quantities
spec.targets.cpasma = init.gdata.cpasma;
spec.weights.cpasma = 1;

spec.targets.li = 1.049;
spec.weights.li = 1;

spec.targets.betap = 0.932;
spec.weights.betap = 1;


% Specify the null locations
spec.targets.rx = [rxP_PRED rxS_PRED];
spec.targets.zx = [zxP_PRED zxS_PRED];

spec.weights.x = [1e6 1e6];

% spec.targets.rbdef = rxPL;
% spec.targets.zbdef = zxPL;
spec.targets.rbdef = rxP_PRED;
spec.targets.zbdef = zxP_PRED;
spec.weights.bdef  = 1;


% Specify the profiles

config.pprime0 = init.gdata.pprime;
config.ffprim0 = init.gdata.ffprim;


% config.pres0 = init.gdata.pres;
% config.fpol0 = init.gdata.fpol;

% pres0,          template pressure profile
% pprime0,        template pprime profile, used if pres0 is missing
% fpol0,          template fpol profile
% ffprim0,        template ffprim profile, used if fpol0 is missing



% Specify coil circuit connections
pp_dir  = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/';
pp_name = 'ppconfig_165286.mat';
load([pp_dir pp_name])
spec.buscode = [0 0 ppconfig.bus_code];

% Load (smoothed) coil currents for the timeslice
load('/u/jwai/d3d_snowflake_2019_wai/coil_currents/165288/coil_currents_smooth.mat')
load('/u/jwai/d3d_snowflake_2019_wai/coil_currents/165288/t.mat')

[~,k] = min(abs(t-time_ms));
cc0 = cc(:,k);

% Define PF coil currents
spec.locks.ic = nan*ones(20,1);

iA = 1:2:17; % index for A coils e.g. F1A, F2A, ..., F9A
iB = 2:2:18; % index for B coils
spec.targets.ic(3:11)  = cc0(iA);
spec.targets.ic(12:20) = cc0(iB);

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


% Run gsdesign
eq = gsdesign(spec, init.gdata, config);

pprime_gs = interp1(eq.psibar, eq.pprime, linspace(0,1,129));
ffprim_gs = interp1(eq.psibar, eq.ffprim, linspace(0,1,129));
jpar_gs   = interp1(eq.psibar, eq.jpar,   linspace(0,1,129));

psibar_gs = linspace(0,1,129);


% Plot the profiles
figure()

subplot(3,1,1)
plot(psibar_efit, pprime_efit, '-b', 'LineWidth', 2)
hold on
plot(psibar_gs, pprime_gs, '-r', 'LineWidth', 2)

grid on
xlabel('Normalized Flux \psi_N')
title('pprime')
    
subplot(3,1,2)
plot(psibar_efit, ffprim_efit, '-b', 'LineWidth', 2)
hold on
plot(psibar_gs, ffprim_gs, '-r', 'LineWidth', 2)

grid on
xlabel('Normalized Flux \psi_N')
title('ffprim')
    
subplot(3,1,3)
plot(psibar_gs, jpar_gs/1e6, '-r', 'LineWidth', 2)

grid on
xlabel('Normalized Flux \psi_N')
title('parallel current density [MA/m^2]')
    
fn = ['eq_cake_constrained_' int2str(shot) '_' int2str(time_ms) '.mat'];
save(fn, 'eq')





