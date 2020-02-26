% function heatsim(eq, shot, time_ms)
shot = 165288;
time_ms = 4000;
load('eq.mat')
if ~exist('eq','var')
    eq = init.gdata;
end

eqdir = '/u/pvail/d3d_snowflake_2019/HeatFluxConstrainedEQ/165288/4000/';
eqnam = 'eq_165288_4000_EQHF.mat';
load([eqdir eqnam])

% ----------------------
% USER / FUNCTION INPUTS
% ----------------------
plotit = 1;
saveit = 1;
saveDir = '/u/jwai/d3d_snowflake_2020/wai/multislice/eq_constrained/';
irConstrained = 1;

% ---------------------------------
% ANALYZE THE SNOWFLAKE EQUILIBRIUM
% ---------------------------------

% Load tokamak definition
tokdir = '/u/jwai/d3d_snowflake_2020/wai/hf_constrained_eq/mat_files/';
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);


% Compute total length of the limiter [m]
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);


% distance to limiter landmarks
s45Deg1 = sLimTot - calcLimDistance(limdata(2,79), limdata(1,79), limdata);
s45Deg2 = sLimTot - calcLimDistance(limdata(2,78), limdata(1,78), limdata);


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
    
rg = eq.rg; 
zg = eq.zg;

% different naming conventions for this info
try   
    bzero = eq.bzero;   % eq generated from gfile
    rzero = eq.rzero;
catch
    bzero = eq.btsurf;  % eq generated from gsdesign
    rzero = eq.rsurf;
end

psizr  = eq.psizr;
psimag = eq.psimag;
psibry = eq.psibry;

% finer mesh
[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);


% Locate the snowflake in the lower divertor
rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);


if abs(psixSL - psibry) < abs(psixPL - psibry)
    swap(psixPL, psixSL);
    swap(rxPL, rxSL);
    swap(zxPL, zxSL);   
end

% determine the snowflake type
snowPlus = false;
snowMinLFS = false;
snowMinHFS = false;

if psixSL > psixPL
    snowPlus = true;
    snowType = 'snowflake plus';
else
    if rxPL < rxSL
        snowMinLFS = true;
        snowType = 'snowflake minus LFS';
    else
        snowMinHFS = true;
        snowType = 'snowflake minus HFS';
    end
end


if plotit    
    contour(rg, zg, psizr, [psixPL psixPL], 'b', 'LineWidth', 2)
    contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)
    
    plot(rxPL, zxPL, 'xb', 'Markersize', 12, 'LineWidth', 3)
    plot(rxSL, zxSL, 'xb', 'Markersize', 12, 'LineWidth', 3)
end


% find the divertor portion of the limiter
Rlessthan = limdata(2,:) <  2.5;
Zlessthan = limdata(1,:) < -0.2;
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

zShelf = -1.25; % DIII-D shelf

% ------------------
% SIMULATE HEAT FLUX
% ------------------

% Model of the midplane SOL

% Determine major radius of primary separatrix on the midplane
C = contourc(rg, zg, psizr, [psixPL psixPL]);
[rmid, k] = max(C(1,:));
zmid = C(2,k);


% Compute the poloidal and total fields at the midplane
[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, rmid, zmid);

BrMid = -1/(2*pi*rmid)*dpsidz;
BzMid =  1/(2*pi*rmid)*dpsidr;

BpMid = sqrt(BrMid*BrMid + BzMid*BzMid);

BtMid = (bzero*rzero)/rmid;

BTotMid = sqrt(BpMid*BpMid + BtMid*BtMid);


% create 1D grid for SOL
lambdaQ = 0.002; % SOL power width [m]
nSOL = 100; % number of grid pts

rSOLMid = linspace(rmid + 5e-5, rmid + 5e-5 + 3*lambdaQ, nSOL)';

psiSOL = interp2(rg, zg, psizr, rSOLMid, zmid);

if plotit
    contour(rg, zg, psizr, [psiSOL(end) psiSOL(end)], '-b')
end

% Compute the power entering the SOL

Pheat = 4; % [MW]
 
fOBL = 0.70; % fraction of power to outboard divertor

PSOL = Pheat;

% power flowing to each divertor leg

PdivOL = PSOL*fOBL;
PdivIL = PSOL*(1-fOBL);

% Compute the peak heat flux at the midplane
q0_parallel_IL = PdivIL/(4*pi*rmid*lambdaQ*(BpMid/BTotMid));    % ADJUST THESE FOR MAGNITUDE ??
q0_parallel_OL = PdivOL/(4*pi*rmid*lambdaQ*(BpMid/BTotMid));


% Compute midplane heat flux profile

dr = rSOLMid - rmid;

q_parallel_IL = (q0_parallel_IL).*exp(-dr/lambdaQ);
q_parallel_OL = (q0_parallel_OL).*exp(-dr/lambdaQ);


%...............................................................................
% Map the midplane heat flux to the divertor

rzLimSOL = zeros(nSOL,4);

% Determine where the flux tubes intersect the limiter

for ii = 1:nSOL
    
    rz = isoflux_spFinder(psizr, psiSOL(ii), rg, zg, limdata, limIdxL);
        
    if snowPlus
        
        idxInner = 1;
        
        switch size(rz,1)
            case 2
                idxOuter = 2;
            case 4
                if any(rz(:,2) == zShelf)
                    idxOuter = find(rz(:,2) == zShelf);
                elseif any(rz(:,1) == limdata(2,74))
                    idxOuter = find(rz(:,1) == limdata(2,74));
                else
                    idxOuter = 4;
                end
            case 6
                if any(rz(:,2) == zShelf)
                    idxOuter = find(rz(:,2) == zShelf);
                end
        end
    
    else
        
        if snowMinLFS
            
            idxInner = 1;
            
            switch size(rz,1)
                case 2
                    idxOuter = 2;
                case 3
                    if rz(2,1) < rxSL
                        idxOuter = 2;
                    else
                        idxOuter = 3;
                    end
                case 4
                    if psiSOL(ii) > psixSL
                        idxOuter = 2;
                    else
                        if any(rz(:,2) == zShelf)
                            idxOuter = find(rz(:,2) == zShelf);
                        elseif rz(2,1) > rxSL
                            idxOuter = 2;
                        else
                            idxOuter = 4;
                        end
                    end
                case 5
                    idxTemp = find(rz(:,1) > rxPL);
                    [~, idxOuter] = max(rz(idxTemp,2));
                    idxOuter = idxTemp(idxOuter);
                case 6
                    [~,idxMax] = max(rz(2:end,2));
                    idxOuter = idxMax + 1;
            end
        
        elseif snowMinHFS
            
            switch size(rz,1)
                case 2
                    idxInner = 1;
                    idxOuter = 2;
                case 3
                    idxInner = 1;
                    idxOuter = 3;
                case 4
                    if psiSOL(ii) > psixSL
                        idxInner = 3;
                        idxOuter = 4;
                    else
                        if any(rz(:,2) == zShelf)
                            idxInner = 1;
                            idxOuter = find(rz(:,2) == zShelf);
                        elseif rz(2,1) > rxPL
                            idxInner = 1;
                            idxOuter = 2;
                        else 
                            idxInner = 1;
                            idxOuter = 4;
                        end
                    end
                case 5
                    idxInner = 1;
                    idxTemp = find(rz(:,1) > rxPL);
                    [~, idxOuter] = max(rz(idxTemp,2));
                    idxOuter = idxTemp(idxOuter);
                case 6
                    idxInner = 1;
                    idxOuter = find(rz(:,2) == zShelf);
                case 7
                    idxInner = 1;
                    idxOuter = find(rz(:,2) == zShelf);
                case 8
                    idxInner = 1;
                    idxOuter = find(rz(:,2) == zShelf);
            end
            
        end
        
    end
    
    rzLimSOL(ii,1) = rz(idxInner,1);
    rzLimSOL(ii,2) = rz(idxInner,2);
    rzLimSOL(ii,3) = rz(idxOuter,1);
    rzLimSOL(ii,4) = rz(idxOuter,2);
    
    if plotit
        plot(rzLimSOL(ii,1), rzLimSOL(ii,2), 'or')
        plot(rzLimSOL(ii,3), rzLimSOL(ii,4), 'or')
    end
            
end



%.....................................................................
% Compute curves which are perpendicular to the flux surfaces (inbd)
[dpsidr, dpsidz] = rzGrad(psizr, rg, zg);
ds = 0.0005;
ns = 5e4;

if snowMinHFS 
    % start the perpendicular curve at the secondary x-pt
    % and then later join this with line between x-pts
    rx = rxSL;
    zx = zxSL;
else
    rx = rxPL;
    zx = zxPL;
end

% seed the initial condition slightly inboard
rperpI(1) = rx - ds;
zperpI(1) = zx + ds;

eta0 = 0.01;
for ii = 1:ns
    
    psi_r = interp2(rg, zg, dpsidr, rperpI(ii), zperpI(ii));
    psi_z = interp2(rg, zg, dpsidz, rperpI(ii), zperpI(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpI(ii+1) = rperpI(ii) - eta*psi_r;
    zperpI(ii+1) = zperpI(ii) - eta*psi_z;
    
    psi = interp2(rg, zg, psizr, rperpI(ii+1), zperpI(ii+1));
    if psi < min(psiSOL)
        break
    end
    
end

rperpI = nonzeros(rperpI);
zperpI = nonzeros(zperpI);
[rerpI, zperpI] = removeNans(rperpI, zperpI);
psiperpI = interp2(rg, zg, psizr, rperpI, zperpI);


%.....................................................................
% Compute curves which are perpendicular to the flux surfaces (outbd)

if snowMinLFS 
    % start the perpendicular curve at the secondary x-pt
    % and then later join this with line between x-pts
    rx = rxSL;
    zx = zxSL;
else
    rx = rxPL;
    zx = zxPL;
end

% seed the initial condition slightly outboard
rperpO(1) = rx + ds;
zperpO(1) = zx - ds;

eta0 = 0.01;
for ii = 1:ns
    
    psi_r = interp2(rg, zg, dpsidr, rperpO(ii), zperpO(ii));
    psi_z = interp2(rg, zg, dpsidz, rperpO(ii), zperpO(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpO(ii+1) = rperpO(ii) - eta*psi_r;
    zperpO(ii+1) = zperpO(ii) - eta*psi_z;
    
    psi = interp2(rg, zg, psizr, rperpO(ii+1), zperpO(ii+1));
    if psi < min(psiSOL)
        break
    end
    
end

rperpO = nonzeros(rperpO);
zperpO = nonzeros(zperpO);
[rerpO, zperpO] = removeNans(rperpO, zperpO);
psiperpO = interp2(rg, zg, psizr, rperpO, zperpO);

if plotit
    plot(rperpI, zperpI, '-b', 'LineWidth', 1)
    plot(rperpO, zperpO, '-b', 'LineWidth', 1)
    
    plot([rxPL rxSL], [zxPL zxSL], '-b', 'LineWidth', 1)
end

 
% interpolate in psi to get rz of divertor entrance
rentrI = interp1(psiperpI, rperpI, psiSOL);
rentrO = interp1(psiperpO, rperpO, psiSOL);
zentrI = interp1(psiperpI, zperpI, psiSOL);
zentrO = interp1(psiperpO, zperpO, psiSOL);

% remove NaNs (produced when interp1 is asked to extrapolate)
[rentrI, rentrO, zentrI, zentrO] = removeNans(rentrI, rentrO, zentrI, ...
    zentrO);

% ^ ^ Do we need to ignore rz entrance that intersect shelf, vert ???


% determine flux tube partitions
idxInner = find(psiSOL > psixSL);  % flux tube MIGHT be between x-points
idxOuter = setdiff(1:nSOL, idxInner); % flux tube outside x-pts
idxShelf = find(rzLimSOL(:,4) == zShelf); % intersect shelf 
idxVert = find(rzLimSOL(:,3) == limdata(2,74));  % intersect outer vert wall
idxBottom = setdiff(1:nSOL, [idxShelf; idxVert]);  % intersect bottom limiter

% find rz coordinates of line between x-points (constant psi-spacing)
rzLine = zeros(length(idxInner),2);

for ii = idxInner
    rzLine(ii,:) = isoflux_cpFinder(psizr, psiSOL(ii), rg, zg,...
        [rxPL rxSL zxPL zxSL]);
end

psiperpX = psiSOL(idxInner);
rentrX = rzLine(:,1);
zentrX = rzLine(:,2);

if plotit
    plot(rentrI, zentrI, 'or') 
    plot(rentrO, zentrO, 'ob') 
    plot(rentrX, zentrX, 'og') 
end

% ................................................
% Compute connection lengths to the divertor target

% inboard
[LdivI,rdivI,zdivI] = calcConnectionLength_wai2(rentrI, zentrI, psizr, rg, zg, ...
    bzero, rzero, limdata, -1);

% outboard
[LdivO,rdivO,zdivO] = calcConnectionLength_wai2(rentrO, zentrO, psizr, rg, zg, ...
    bzero, rzero, limdata, 1);
    
% line between xpoints
if ~snowPlus
    dir = 1;       % outbd
    if snowMinHFS  % inbd
        dir = -1;
    end
    
    [LdivX, rdivX, zdivX] = calcConnectionLength_wai2(rentrX, zentrX, ...
        psizr, rg, zg, bzero, rzero, limdata, dir);
end

if plotit
    plot(rdivI, zdivI, 'or')
    plot(rdivO, zdivO, 'ob')
    if ~snowPlus, plot(rdivX, zdivX, 'og'); end
end


%.........................
% Heat equation parameters

chi = 0.01; % thermal  diffusivity [Wb^2/s]
Ti = 50 * 11600;  % [K]   Ion and electron temperatures
Te = 50 * 11600;  % [K]
k = 1.38e-23;  % Boltzmann constant [(m^2*kg)/(s^2*K)]
mi = 3.344e-27; % Deuterium mass [kg]
cs = sqrt(k*(Ti + Te)/mi);  % Sound speed

% Compute the time-of-flight for each HFS/LFS fieldline
tauI = LdivI/cs;
tauO = LdivO/cs;
if ~snowPlus, tauX = LdivX/cs; end



%........................................................
% Compute the heat flux profile region between 2 x-points
if ~snowPlus
    nRegion = length(LdivX);
    iRegion = 1:nRegion;

    frad = 0.80; 
    
    if snowMinHFS
        q_par_midplane = q_parallel_IL(iRegion); 
    elseif snowMinLFS
        q_par_midplane = q_parallel_OL(iRegion);
    end
    
    psiSOLX = psiSOL(iRegion);
    
    % find heat flux
    [qdiv_perpX, sDivX, qdiv_parX] = ...
        find_qperp(nRegion, iRegion, psiperpX, tauX, q_par_midplane, frad, ...
        chi, psiSOLX, rdivX, zdivX, limdata, psizr, rg, zg, bzero, rzero, ...
        sLimTot);
end


%........................................................
% Compute the heat flux profile for the inboard/HFS region

% inboard parameters
nRegion = length(LdivI);

if  snowPlus || snowMinLFS
    iRegion = 1:nRegion;
elseif snowMinHFS
    iRegion = length(LdivX):length(LdivX) + nRegion;    
end

q_par_midplane = q_parallel_IL(iRegion);
psiSOLI = psiSOL(iRegion);
frad = 0.80; 

% find heat flux
[qdiv_perpI, sDivI, qdiv_parI] = ...
    find_qperp(nRegion, iRegion, psiperpI, tauI, q_par_midplane, frad, ...
    chi, psiSOLI, rdivI, zdivI, limdata, psizr, rg, zg, bzero, rzero, ...
    sLimTot);


%........................................................
% Compute the heat flux profile for the outboard/LFS region
nRegion = length(LdivO);

if snowPlus || snowMinHFS
    iRegion = 1:nRegion;
elseif snowMinLFS
    iRegion = length(LdivX)+1:length(LdivX) + nRegion;
end

q_par_midplane = q_parallel_OL(iRegion);
psiSOLO = psiSOL(iRegion);
frad = 0.80; 

% find heat flux
[qdiv_perpO, sDivO, qdiv_parO] = ...
    find_qperp(nRegion, iRegion, psiperpO, tauO, q_par_midplane, frad, ...
    chi, psiSOLO, rdivO, zdivO, limdata, psizr, rg, zg, bzero, rzero, ...
    sLimTot);


%............................
% Find strike point locations 

[~,iIP] = min(abs(psiSOLI - psixPL));  % index of pt in the inboard region 
                                       % closest to psi of the primary x-pt                                       
[~,iIS] = min(abs(psiSOLI - psixSL));  % index of pt " " " secondary x-pt
[~,iXP] = min(abs(psiSOLX - psixPL));
[~,iXS] = min(abs(psiSOLX - psixSL));
[~,iOP] = min(abs(psiSOLO - psixPL));
[~,iOS] = min(abs(psiSOLO - psixSL));

if snowPlus
    iI = iIP;
    iO = iOP;
    iX = NaN;
elseif snowMinLFS
    iI = iIP;
    iO = iOS;
    iX = iXP;
elseif snowMinHFS
    iI = iIS;
    iO = iOP;
    iX = iXP;
end

[sSPI, sSPX, sSPO] = deal(sDivI(iI), sDivX(iX), sDivO(iO));
[spRI, spRX, spRO] = deal(rdivI(iI), rdivX(iX), rdivO(iO));
[spZI, spZX, spZO] = deal(zdivI(iI), zdivX(iX), zdivO(iO));


% ..............
% plot heat flux
if plotit
    
   figure(2)
   hold on
   
   xlabel('s [cm]')
   ylabel('Heat Flux [MW/m^2]')
   title([int2str(shot) ': ' int2str(time_ms) ' ms'])

   maxQ = max([qdiv_perpI; qdiv_perpX; qdiv_perpO]);   
   ylim(1.10*[0 maxQ])
   
   xline(sSPI, 'k');
   xline(sSPO, 'k');
   xline(sSPX, 'k');
   xline(s45Deg1, '--k');
   xline(s45Deg2, '--k');
   
   plot(sDivX, qdiv_perpX, '-og', 'LineWidth', 1, 'MarkerSize', 2)
   plot(sDivI, qdiv_perpI, '-or', 'LineWidth', 1, 'MarkerSize', 2)
   plot(sDivO, qdiv_perpO, '-ob', 'LineWidth', 1, 'MarkerSize', 2)   
end

%...............................................................................
% Fit the Eich profile to the simulated heat flux data

qEich_Inner = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (-(x-s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x-s0))./S) + qBG;

qEich_Outer = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;

ftInner = fittype(qEich_Inner, 'problem', 'fExp', 'independent', 'x');
ftOuter = fittype(qEich_Outer, 'problem', 'fExp', 'independent', 'x');

optionsInner = fitoptions(ftInner);
optionsInner.Lower = [0 0 0 0 0];

optionsOuter = fitoptions(ftOuter);
optionsOuter.Lower = [0 0 0 0 0];


%............................
% Eich profile fit to INBOARD

fExp = calcFluxExpansion(spRI, spZI, psizr, rg, zg, psixPL, bzero, rzero);

optionsInner.StartPoint = [2*max(qdiv_perpI) 0.01 lambdaQ sSPI 0];
   
fitEichI = fit(sDivI, qdiv_perpI, ftInner, 'problem', fExp, optionsInner);

if plotit    
    s_Eich = linspace(1,1.7,1000);
    
    qperp_EichI = eichProfileInner(s_Eich, fitEichI.q0, ...
        fitEichI.S, fitEichI.lambdaQ, fExp, fitEichI.s0);
    
    figure(2)
    plot(s_Eich, qperp_EichI, '-k', 'LineWidth', 1)
end


%............................
% Eich profile fit to X-PTS REGION
if ~snowPlus    
    fExp = calcFluxExpansion(spRX, spZX, psizr, rg, zg, psixPL, bzero, rzero);

    if snowMinLFS                
    
        optionsOuter.StartPoint = [2*max(qdiv_perpX) 0.01 lambdaQ sSPX 0];
        
        fitEichX = fit(sDivX, qdiv_perpX, ftOuter, 'problem', fExp, optionsOuter);            
        
        qperp_EichX = eichProfileOuter(s_Eich, fitEichX.q0, ...
            fitEichX.S, fitEichX.lambdaQ, fExp, fitEichX.s0);        

    elseif snowMinHFS        
        
        optionsInner.StartPoint = [2*max(qdiv_perpX) 0.01 lambdaQ sSPX 0];
        
        fitEichX = fit(sDivX, qdiv_perpX, ftInner, 'problem', fExp, optionsInner);
        
        qperp_EichX = eichProfileInner(s_Eich, fitEichX.q0, ...
            fitEichX.S, fitEichX.lambdaQ, fExp, fitEichX.s0);        
    end
                       
    if plotit
        s_Eich = linspace(1,1.7,1000);                
        figure(2)
        plot(s_Eich, qperp_EichX, '-k', 'LineWidth', 1)
    end
end


%............................
% Eich profile fit to OUTBOARD

fExp = calcFluxExpansion(spRO, spZO, psizr, rg, zg, psixPL, bzero, rzero);

optionsOuter.StartPoint = [2*max(qdiv_perpO) 0.01 lambdaQ sSPO 0];
   
fitEichO = fit(sDivO, qdiv_perpO, ftOuter, 'problem', fExp, optionsOuter);

if plotit    
    s_Eich = linspace(1,1.7,1000);
    
    qperp_EichO = eichProfileOuter(s_Eich, fitEichO.q0, ...
        fitEichO.S, fitEichO.lambdaQ, fExp, fitEichO.s0);
    
    figure(2)
    plot(s_Eich, qperp_EichO, '-k', 'LineWidth', 1)
end


%...............................................................................
% Save the data   
hfsim = struct(              ...
    'shot',          shot,           ...
    'time_ms',       time_ms,        ...
    'snowType',      snowType,       ...
    'rg257',         rg,             ...
    'zg257',         zg,             ...
    'psizr257257',   psizr,          ...
    'rxPL',          rxPL,           ...
    'rxSL',          rxSL,           ...
    'zxPL',          zxPL,           ...
    'zxSL',          zxSL,           ...
    'psixPL',        psixPL,         ...
    'psixSL',        psixSL,         ...
    'rentrI',        rentrI,         ...
    'zentrI',        zentrI,         ...
    'rentrX',        rentrX,         ...
    'zentrX',        zentrX,         ...
    'rentrO',        rentrO,         ...
    'zentrO',        zentrO,         ...
    'LdivI',         LdivI,          ...
    'LdivX',         LdivX,          ...
    'LdivO',         LdivO,          ...    
    'rdivI',         rdivI,          ...
    'zdivI',         zdivI,          ...  
    'rdivX',         rdivX,          ...
    'zdivX',         zdivX,          ...  
    'rdivO',         rdivO,          ...
    'zdivO',         zdivO,          ...      
    'spRI',          spRI,           ...
    'spRX',          spRX,           ... 
    'spRO',          spRO,           ...
    'spZI',          spZI,           ...
    'spZX',          spZX,           ... 
    'spZO',          spZO,           ...
    'sSPI',          sSPI,           ...
    'sSPX',          sSPX,           ...
    'sSPO',          sSPO,           ...
    'qdiv_parI',     qdiv_parI,      ...
    'qdiv_parX',     qdiv_parX,      ...
    'qdiv_parO',     qdiv_parO,      ...
    'qdiv_perpI',    qdiv_perpI',    ...
    'qdiv_perpX',    qdiv_perpX',    ...
    'qdiv_perpO',    qdiv_perpO',    ... 
    'fitEichI',      fitEichI,       ...
    'fitEichX',      fitEichX,       ...
    'fitEichO',      fitEichO,       ...
    's45Deg1',       s45Deg1,        ...
    's45Deg2',       s45Deg2         ...
);

if saveit
    
    % save data
    if irConstrained        
        fn = ['hfsim_c_' int2str(shot) '_' int2str(time_ms) '.mat']; 
    else
        fn = ['hfsim_u_' int2str(shot) '_' int2str(time_ms) '.mat']; 
    end        
    save([saveDir fn], 'hfsim')
    
    % save figs
    fn = ['sfd_geo_' int2str(shot) '_' int2str(time_ms)]; 
    figure(11)
    savefig([saveDir fn]);
    
    fn = ['heat_profile_' int2str(shot) '_' int2str(time_ms)]; 
    figure(2)
    savefig([saveDir fn]);       
end












