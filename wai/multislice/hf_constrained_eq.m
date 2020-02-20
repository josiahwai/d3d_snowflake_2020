% Functionalized versions of Pat's code. 
% 
% - obtain heat flux parameters
% - predict xpoint locations
% - design new equilibrium
% 
% - measure edge current
% - simulate heat flux and compare


% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
shot = 165288;
time_ms = 4000;  

plotit = 1;

saveit = 0;

simulate_heat = 1;


% ------------------------------------------
% OBTAIN SNOWFLAKE AND HEAT FLUX PARAMETERS
% ------------------------------------------

addpath(genpath('/u/jwai/d3d_snowflake_2019_wai'));
addpath(genpath('/u/jwai/d3d_snowflake_2019'));


% Load tokamak definition
load('/u/jwai/d3d_snowflake_2019_wai/hf_constrained_eq/mat_files/d3d_obj_mks_struct_129129.mat')


% Compute total length of the limiter [m]
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% Configure the plots
if plotit    
    figure(11)
    plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
    hold on
    axis equal
    axis([1.0 1.5 -1.4 -0.9])
    xlabel('R [m]')
    ylabel('Z [m]')
    title([int2str(shot) ': ' int2str(time_ms) ' ms'])
end
    

% Load the equilibrium
eq_dir = strcat('/u/jwai/d3d_snowflake_2019_wai/multislice/eq_unconstrained/', int2str(shot), '/cake');
eq = read_eq(shot, time_ms/1000, eq_dir);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

rmaxis = eq.gdata.rmaxis;
zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
psimag = eq.gdata.psimag;
psibry = eq.gdata.psibry;

[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);


% Locate the snowflake in the lower divertor

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    
    temp   = psixPL;
    psixPL = psixSL;
    psixSL = temp;
    
    temp1 = rxPL;
    temp2 = zxPL;
    rxPL  = rxSL;
    zxPL  = zxSL;
    rxSL  = temp1;
    zxSL  = temp2;
    
end

if plotit    
    contour(rg, zg, psizr, [psixPL psixPL], 'b', 'LineWidth', 2)
    contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)
    
    plot(rxPL, zxPL, 'xb', 'Markersize', 12, 'LineWidth', 3)
    plot(rxSL, zxSL, 'xb', 'Markersize', 12, 'LineWidth', 3)
end

% Find the strike points in the EFIT

Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);

r_SP1_EFIT01 = spRZLP(1,1); 
r_SP2_EFIT01 = spRZLP(2,1); 
r_SP3_EFIT01 = spRZLS(4,1);

z_SP1_EFIT01 = spRZLP(1,2); 
z_SP2_EFIT01 = spRZLP(2,2); 
z_SP3_EFIT01 = spRZLS(4,2);

if plotit    
    plot(r_SP1_EFIT01, z_SP1_EFIT01, 'ob', 'LineWidth', 2)
    plot(r_SP2_EFIT01, z_SP2_EFIT01, 'ob', 'LineWidth', 2)
    plot(r_SP3_EFIT01, z_SP3_EFIT01, 'ob', 'LineWidth', 2)
end

 
fExp_SP1_EFIT01 = calcFluxExpansion(r_SP1_EFIT01, z_SP1_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);
fExp_SP2_EFIT01 = calcFluxExpansion(r_SP2_EFIT01, z_SP2_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);
fExp_SP3_EFIT01 = calcFluxExpansion(r_SP3_EFIT01, z_SP3_EFIT01, psizr, ...
    rg, zg, psibry, bzero, rzero);

s_SP1_EFIT01 = sLimTot - calcLimDistance(r_SP1_EFIT01, z_SP1_EFIT01, limdata); 
s_SP2_EFIT01 = sLimTot - calcLimDistance(r_SP2_EFIT01, z_SP2_EFIT01, limdata);
s_SP3_EFIT01 = sLimTot - calcLimDistance(r_SP3_EFIT01, z_SP3_EFIT01, limdata);


%....................................
% Load the heat flux data for the shot

qperp_dir  = ['/u/jwai/d3d_snowflake_2019_wai/multislice/ir_data/'];
qperp_data = [num2str(shot) '/qperp_' num2str(shot) '.mat'];
sFile = 's.mat';

load([qperp_dir qperp_data])
load([qperp_dir sFile]) 

[~,k] = min(abs(t-time_ms));
qperp = qperp_all(k,:)';


% Remove the gap

idxGap = find(s < 170);

% Remove the gap

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

if plotit 
    figure(22)
    plot(s, qperp, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
    hold on
    plot([50 210], [0 0], '--k')
    
    axis([80 180 -2 50])
    xlabel('s [cm]')
    ylabel('Heat Flux W/cm^2')
    title([int2str(shot) ': ' int2str(time_ms) ' ms'])
end


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
rtree_dir = '/u/jwai/d3d_snowflake_2019/invMap/tree/v6/';
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

config = tok_data_struct;
config.max_iterations = 30;

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

spec.targets.rsep = init.gdata.rbbbs([1:60 85:89]);
spec.targets.zsep = init.gdata.zbbbs([1:60 85:89]);
spec.weights.sep = 1e3*ones(length(spec.targets.rsep),1);

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
spec.targets.bdef  = 1;

%.....................
% Specify the profiles

config.constraints = 1;  % allow for scaling/peaking of profiles
config.pres0 = init.gdata.pres;
config.fpol0 = init.gdata.fpol;

%.................................
% Specify coil circuit connections
pp_dir = ['/u/jwai/d3d_snowflake_2019_wai/multislice/coil_data/' num2str(shot)];
pp_data = ['/ppconfig_' num2str(shot) '.mat'];

load([pp_dir pp_data])
spec.buscode = [0 0 ppconfig.bus_code];


%........................
% Load coil currents
coil_dir = ['/u/jwai/d3d_snowflake_2019_wai/multislice/coil_data/' num2str(shot)];
coil_data = '/coil_currents_smooth.mat';
t_data = '/t.mat';

load([coil_dir coil_data]);
load([coil_dir t_data]);

[~,k] = min(abs(t-time_ms));
cc0 = cc(:,k);


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
% Run gsdesign
spec.fig = -1;

eq = gsdesign(spec, init.gdata, config);

out_dir = '/u/jwai/d3d_snowflake_2019_wai/multislice/eq_constrained/';
out_name = ['eq_' num2str(shot) '_' num2str(time_ms) '_constrained'];

save([out_dir out_name], 'eq');



% --------------------
% COMPARE EDGE CURRENT
% --------------------

jpar = eq.jpar;
psin = eq.psibar;


% load previous edge current
jdir = '/u/jwai/d3d_snowflake_2019_wai/multislice/eq_unconstrained/165288/';
jdata = 'jdata.mat';

jpar0 = jtot;
psin0 = psi_n;



if plotit
    figure(33)
    hold on
    plot(psin,  jpar,  'b', 'linewidth', 2)
    plot(psin0, jpar0, 'r', 'linewidth', 2)
end


% --------------------
% SIMULATE HEAT FLUX
% --------------------

if simulate_heat
    heatsim
end





























