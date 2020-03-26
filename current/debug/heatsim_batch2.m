function heatsim_batch2(eq, shot, time_ms, opts)

% ---------------------------------
% ANALYZE THE SNOWFLAKE EQUILIBRIUM
% ---------------------------------
close all
struct_to_ws(opts);

% Load tokamak definition
tokdir = [root 'inputs/tok_data/'];
tokdata = 'd3d_obj_mks_struct_129129.mat';
load([tokdir tokdata]);


% Compute total length of the limiter [m]
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% interpolate for higher point density along limiter 
[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 500, true);

% distance to limiter landmarks
s45Deg1 = sLimTot - calcLimDistance(limdata(2,79), limdata(1,79), limdata);
s45Deg2 = sLimTot - calcLimDistance(limdata(2,78), limdata(1,78), limdata);


% Configure the plots
if plotit    
    figure(11)
    plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
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
    % check the sign of the angle defined by the secondary x-pt, primary x-pt, and
    % plasma axis (with primary x-pt in the middle) to determine HFS/LFS
    % This is more robust than comparing rxPL, rxSL
    
    a = [eq.rmaxis eq.zmaxis 0] - [rxPL zxPL 0];
    b = [rxSL zxSL 0] - [rxPL zxPL 0];
    
    c = cross(a,b);
    
    if c(3) < 0  % CW rotation from a to b
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


%...........................................
% Map the midplane heat flux to the divertor


%...................................................................
% Compute curves which are perpendicular to the flux surfaces (inbd)
[dpsidr, dpsidz] = rzGrad(psizr, rg, zg);
ds = .005;
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
[rperpI, zperpI] = removeNans(rperpI, zperpI);
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
[rperpO, zperpO] = removeNans(rperpO, zperpO);
psiperpO = interp2(rg, zg, psizr, rperpO, zperpO);

if plotit
    plot(rperpI, zperpI, '-b', 'LineWidth', 1)
    plot(rperpO, zperpO, '-b', 'LineWidth', 1)
    
    plot([rxPL rxSL], [zxPL zxSL], '-b', 'LineWidth', 1)
end

 
% interpolate in psi to get psi-spaced rz at divertor entrance
rentrI = interp1(psiperpI, rperpI, psiSOL);
rentrO = interp1(psiperpO, rperpO, psiSOL);
zentrI = interp1(psiperpI, zperpI, psiSOL);
zentrO = interp1(psiperpO, zperpO, psiSOL);

% remove NaNs (produced when interp1 is asked to extrapolate)
[rentrI, rentrO, zentrI, zentrO] = removeNans(rentrI, rentrO, zentrI, ...
    zentrO);

% TAKE WARNING ^ ^: may need to ignore rz entrance that intersect shelf, vert

% determine flux tube partitions
idxInner = find(psiSOL > psixSL);  % flux tube MIGHT be between x-points
idxOuter = setdiff(1:nSOL, idxInner); % flux tube outside x-pts


% find rz coordinates of line between x-points (constant psi-spacing)
rzLine = zeros(length(idxInner),2);

for ii = idxInner
    rzLine(ii,:) = isoflux_cpFinder(psizr, psiSOL(ii), rg, zg,...
        [rxPL rxSL zxPL zxSL]);
end

psiperpX = psiSOL(idxInner);
rentrX = rzLine(:,1);
zentrX = rzLine(:,2);

% x-points are so close the psi-spacing is negligible ==> perfect snowflake
if length(rentrX) < 1
    perfectSnow = true;  
else
    perfectSnow = false;
end

if plotit
    plot(rentrI, zentrI, 'or') 
    plot(rentrO, zentrO, 'ob') 
    plot(rentrX, zentrX, 'og') 
end


% ................................................
% Compute connection lengths to the divertor target

% Positions at divertor entrance (rentr, zentr) as found may lie outside
% the limiter. Check if this is the case and if so, project the divertor 
% entrance back to the limiter by integrating connection length backwards. 


% inboard
iOut = ~inpolygon(rentrI,zentrI,rlim,zlim);
startinI(1:length(rentrI)) = true;
startinI(iOut) = false;

[LdivI,rdivI,zdivI] = calcConnectionLength(rentrI, zentrI, psizr, rg, zg, ...
    bzero, rzero, limdata, -1, startinI);

LdivI(iOut) = 0;

% outboard
iOut = ~inpolygon(rentrO,zentrO,rlim,zlim);
startinO(1:length(rentrO)) = true;
startinO(iOut) = false;

[LdivO,rdivO,zdivO] = calcConnectionLength(rentrO, zentrO, psizr, rg, zg, ...
    bzero, rzero, limdata, 1, startinO);

LdivO(iOut) = 0;

    
% line between xpoints
LdivX = NaN; 
rdivX = NaN; 
zdivX = NaN;
if ~snowPlus && ~perfectSnow
    dir = 1;       % outbd
    if snowMinHFS  % inbd
        dir = -1;
    end
    
    iOut = ~inpolygon(rentrX,zentrX,rlim,zlim);
    startinX(1:length(rentrX)) = true;
    startinX(iOut) = false;
    
    [LdivX, rdivX, zdivX] = calcConnectionLength(rentrX, zentrX, ...
        psizr, rg, zg, bzero, rzero, limdata, dir, startinX);
    
    LdivX(iOut) = 0;
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
qdiv_perpX = NaN;
sDivX = NaN;
psiSOLX = NaN;
qdiv_parX = NaN;
if ~snowPlus && ~perfectSnow
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
    iRegion = length(LdivX)+1:length(LdivX) + nRegion;    
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
                                       % closest to psi at primary x-pt                                       
[~,iIS] = min(abs(psiSOLI - psixSL));  % index of pt " " " secondary x-pt
[~,iXP] = min(abs(psiSOLX - psixPL));
[~,iOP] = min(abs(psiSOLO - psixPL));
[~,iOS] = min(abs(psiSOLO - psixSL));

if snowPlus
    iI = iIP;
    iO = iOP;
    iX = 1;
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
   xline(s45Deg1, '--k');
   xline(s45Deg2, '--k');
   
   plot(sDivI, qdiv_perpI, '-or', 'LineWidth', 1, 'MarkerSize', 2)
   plot(sDivO, qdiv_perpO, '-ob', 'LineWidth', 1, 'MarkerSize', 2)   

   if ~snowPlus && ~perfectSnow
       xline(sSPX, 'k');
       plot(sDivX, qdiv_perpX, '-og', 'LineWidth', 1, 'MarkerSize', 2)
   end
   
   % .............
   % IR-based data
   
   % Load heat flux data q(s,t), s=distance along limiter, and t=time
   qperp_dir  = [root 'inputs/qperp/' num2str(shot) '/'];
   qperp_data = ['qperp_' num2str(shot) '.mat'];
   
   load([qperp_dir qperp_data])  % loads q, s, and t
   
   [~,k] = min(abs(t-time_ms));
   qperp = qperp(k,:)'/500;
   s = s/100;
   
   % Remove the gap from s (distance along limiter)
   iGap = find(s < 1.70,1,'last');
   dgap = s(iGap+1) - s(iGap);
   
   s(iGap(end)+1:end) = s(iGap(end)+1:end) - dgap;
   
   plot(s, qperp, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
   
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
    's45Deg1',       s45Deg1,        ...
    's45Deg2',       s45Deg2         ...
);

if saveit
    fn = ['hfsim_' int2str(shot) '_' int2str(time_ms) '.mat']; 
    save([saveDir fn], 'hfsim')
    
    % save figs
    fn = ['sfd_geo_' int2str(shot) '_' int2str(time_ms)]; 
    figure(11); h = gcf;
    savefig(h, [saveDir fn]);
    
    fn = ['heat_profile_' int2str(shot) '_' int2str(time_ms)]; 
    figure(2); h = gcf;
    savefig(h, [saveDir fn]);       
end















