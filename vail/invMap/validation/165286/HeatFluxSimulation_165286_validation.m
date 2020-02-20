
time = 2000;

shot = 165286;

plotit = 1;
saveit = 0;
 
% shotdir = ['shotdir/' int2str(shot) '/efit03/'];

% geqdsk_name = ['g' int2str(shot) '.0' int2str(time)];

%...............................................................................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

limdata = limdata(:,[1:79 81:end]);

% Compute total length of the limiter [m]

sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

% Compute distance to limiter landmarks [m]

s45Deg1 = sLimTot - calcLimDistance(limdata(2,79), limdata(1,79), limdata);
s45Deg2 = sLimTot - calcLimDistance(limdata(2,78), limdata(1,78), limdata);

%...............................................................................
% Configure the plots

if plotit
    plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
    hold on
    axis equal
    axis([1.0 1.5 -1.4 -0.9])
    xlabel('R [m]')
    ylabel('Z [m]')
    title([int2str(shot) ': ' int2str(time) ' ms'])
end

%...............................................................................
% Load and analyze the standard EFIT

efit_dirname = '/u/pvail/d3d_snowflake_2019/invMap/validation/165286';
    
eq = read_eq(shot, time/1000, efit_dirname);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
psibry = eq.gdata.psibry;

[psizr257257, rg257, zg257] = regrid(rg, zg, psizr, 257, 257);

% Compute Br and Bz at the grid points

[dpsidr, dpsidz] = rzGrad(psizr257257, rg257, zg257);

Br257257 = zeros(size(psizr257257,1),size(psizr257257,2));
Bz257257 = zeros(size(psizr257257,1),size(psizr257257,2));

for ii = 1:length(rg257)
    for jj = 1:length(zg257)
        
        Br257257(jj,ii) = -1/(2*pi*rg257(ii))*dpsidz(jj,ii);
        Bz257257(jj,ii) =  1/(2*pi*rg257(ii))*dpsidr(jj,ii);
        
    end
end

% Analyze the equilibrium

% Lower divertor (SFD)

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr257257, rExp, zExp, rhoExp, rg257, zg257);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr257257, rxPL, zxPL, rg257, zg257);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr257257, rxSL, zxSL, rg257, zg257);

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
    contour(rg257, zg257, psizr257257, [psixPL psixPL], 'k', 'LineWidth', 2)
    contour(rg257, zg257, psizr257257, [psixSL psixSL], 'k', 'LineWidth', 1)
    
    plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)
    plot(rxSL, zxSL, 'xk', 'Markersize', 12, 'LineWidth', 3)
end

% Find the strike points in the standard EFIT

% Define strike point segments
 
Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

zShelf = limdata(1,limIdxL(3));

% Lower divertor (SFD)

spRZLP = isoflux_spFinder(psizr257257, psixPL, rg257, zg257, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr257257, psixSL, rg257, zg257, limdata, limIdxL);

% Determine which SPs are the true strike points

if psixSL > psixPL % SFD-Plus
    
    if size(spRZLP,1) == 3
        spRZLP = spRZLP([1 3],:);
    elseif size(spRZLP,1) == 4
        spRZLP = spRZLP([1 4],:);
    end
    if size(spRZLS,1) == 4
        spRZLS = spRZLS([2 3],:);
    end
    
else % SFD-Minus
    
end

if plotit
    for ii = 1:size(spRZLP,1)
        plot(spRZLP(ii,1), spRZLP(ii,2), 'ob', 'LineWidth', 2)
    end
    for ii = 1:size(spRZLS,1)
        plot(spRZLS(ii,1), spRZLS(ii,2), 'om', 'LineWidth', 2)
    end
end

% Compute distance to strike points [m]

sSPP1 = sLimTot - calcLimDistance(spRZLP(1,1), spRZLP(1,2), limdata); 
sSPP2 = sLimTot - calcLimDistance(spRZLP(2,1), spRZLP(2,2), limdata);

sSPS1 = sLimTot - calcLimDistance(spRZLS(1,1), spRZLS(1,2), limdata); 
sSPS2 = sLimTot - calcLimDistance(spRZLS(2,1), spRZLS(2,2), limdata);

%...............................................................................
% Modeling of the midplane SOL

% Determine major radius of primary and secondary separatrix on the midplane
 
psimid = interp2(rg257, zg257, psizr257257, rg257, zmaxis);
 
idxP = find(psimid > psixPL, 1, 'last' ); % midplane primary
idxS = find(psimid > psixSL, 1, 'last' ); % midplane secondary
 
idxPrimary   = [idxP-1 idxP idxP+1 idxP+2];
idxSecondary = [idxS-1 idxS idxS+1 idxS+2];

ppPrimary   = spline(rg257(idxPrimary),   psimid(idxPrimary));
ppSecondary = spline(rg257(idxSecondary), psimid(idxSecondary));

cP = ppPrimary.coefs(2,:);
cS = ppSecondary.coefs(2,:);

rootsP = roots([cP(1) cP(2) cP(3) cP(4)-psixPL]) + rg257(idxP);
rootsS = roots([cS(1) cS(2) cS(3) cS(4)-psixSL]) + rg257(idxS);

[~,idxminP] = min(abs(rootsP-rg257(idxP)));
[~,idxminS] = min(abs(rootsS-rg257(idxS)));

rmidPrimary   = rootsP(idxminP);
rmidSecondary = rootsS(idxminS);

% Radial grid at the plasma edge

lambdaQ = 0.002; % SOL power width [m]
 
rSOLMid = linspace(rmidPrimary, rmidPrimary + 3*lambdaQ, 100)';

% Compute flux at each grid point

psiSOL = interp2(rg257, zg257, psizr257257, rSOLMid, zmaxis);

if plotit
    contour(rg257, zg257, psizr257257, [psiSOL(end) psiSOL(end)], '-b')
end

% Compute the ratio of poloidal to total magnetic field at the midplane

BrMid = interp2(rg257, zg257, Br257257, rmidPrimary, zmaxis);
BzMid = interp2(rg257, zg257, Bz257257, rmidPrimary, zmaxis);

BpMid = sqrt(BrMid*BrMid + BzMid*BzMid);

BTMid = (bzero*rzero)/rmidPrimary;

BTotMid = sqrt(BrMid*BrMid + BzMid*BzMid + BTMid*BTMid);

% Compute the power entering the SOL

Pheat = 4; % [MW]

frad = 0.80; % radiated power fraction
fOBL = 0.55; % fraction of power to outboard divertor

PSOL = Pheat*(1-frad);

% Compute the power flowing to each divertor leg

PdivOL = PSOL*fOBL;
PdivIL = PSOL*(1-fOBL);

% Compute the peak heat flux at the midplane

q0_parallel_IL = PdivIL/(4*pi*rmidPrimary*lambdaQ*(BpMid/BTotMid));
q0_parallel_OL = PdivOL/(4*pi*rmidPrimary*lambdaQ*(BpMid/BTotMid));

% Compute midplane heat flux profile

dr = rSOLMid - rmidPrimary;

q_parallel_IL = (q0_parallel_IL).*exp(-dr/lambdaQ);
q_parallel_OL = (q0_parallel_OL).*exp(-dr/lambdaQ);

%...............................................................................
% Map the midplane heat flux to the divertor

rzLimSOL = zeros(length(psiSOL),4);

% Determine where the flux tubes intersect the limiter

for ii = 1:length(psiSOL)
    
    rz = isoflux_spFinder(psizr257257, psiSOL(ii), rg257, zg257, ...
        limdata, limIdxL);
        
    if psixSL > psixPL % SFD-Plus
        
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
    
    else % SFD-Minus
         
    end
    
    rzLimSOL(ii,1) = rz(idxInner,1);
    rzLimSOL(ii,2) = rz(idxInner,2);
    rzLimSOL(ii,3) = rz(idxOuter,1);
    rzLimSOL(ii,4) = rz(idxOuter,2);
               
end

% Determine the flux tube partitioning

idxInner = find(psiSOL > psixSL);
%idxOuter = setdiff(1:length(psiSOL), idxInner);

% Determine the indices of flux tubes which intersect the shelf

idxShelf = find(rzLimSOL(:,4) == zShelf);

% Determine the indices of flux tubes which intersect outer vertical wall

idxVert = find(rzLimSOL(:,3) == limdata(2,74));

% Determine the indices which intersect the bottom limiter

idxBottom = setdiff(1:length(psiSOL), [idxShelf; idxVert]);

%...............................................................................
% Compute the connection lengths for SOL fieldlines

% Compute curves which are perpendicular to the flux surfaces (primary)

ns = 5e4;

rperpP = zeros(ns+1,1);
zperpP = zeros(ns+1,1);

rperpP(1) = 1.001*rxPL;
zperpP(1) = 1.001*zxPL;

eta0 = 0.01;

for ii = 1:ns
    
    psi_r = interp2(rg257, zg257, dpsidr, rperpP(ii), zperpP(ii));
    psi_z = interp2(rg257, zg257, dpsidz, rperpP(ii), zperpP(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpP(ii+1) = rperpP(ii) - eta*psi_r;
    zperpP(ii+1) = zperpP(ii) - eta*psi_z;
    
    psi = interp2(rg257, zg257, psizr257257, rperpP(ii+1), zperpP(ii+1));
    if abs((psi - psiSOL(end))) < 1e-5
        break
    end
    
end

rperpP = nonzeros(rperpP);
zperpP = nonzeros(zperpP);

psiperpP = interp2(rg257, zg257, psizr257257, rperpP, zperpP);

% Compute curves which are perpendicular to the flux surfaces (secondary)

ns = 5e4;

rperpS = zeros(ns+1,1);
zperpS = zeros(ns+1,1);

rperpS(1) = 0.999*rxPL;
zperpS(1) = 1.001*zxPL;

eta0 = 0.01;

for ii = 1:ns
    
    psi_r = interp2(rg257, zg257, dpsidr, rperpS(ii), zperpS(ii));
    psi_z = interp2(rg257, zg257, dpsidz, rperpS(ii), zperpS(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpS(ii+1) = rperpS(ii) - eta*psi_r;
    zperpS(ii+1) = zperpS(ii) - eta*psi_z;
    
    psi = interp2(rg257, zg257, psizr257257, rperpS(ii+1), zperpS(ii+1));
    if abs((psi - psiSOL(end))) < 1e-5
        break
    end
    
end

rperpS = nonzeros(rperpS);
zperpS = nonzeros(zperpS);

psiperpS = interp2(rg257, zg257, psizr257257, rperpS, zperpS);

if plotit
    plot(rperpP, zperpP, '-b', 'LineWidth', 1)
    plot(rperpS, zperpS, '-b', 'LineWidth', 1)
    
end

% Compute intersection points of flux tubes with the divertor entrance
 
rentr = zeros(length(psiSOL),2);
zentr = zeros(length(psiSOL),2);
 
for ii = 1:length(psiSOL)     
    
    psiSOLi = psiSOL(ii);
    
    if psixSL > psixPL % SFD-Plus
        
        % HFS SOL 
        
        [~,idx] = min(abs(psiSOLi - psiperpS));
        
        rentr(ii,1) = rperpS(idx);
        zentr(ii,1) = zperpS(idx);
        
        % LFS SOL
        
        [~,idx] = min(abs(psiSOLi - psiperpP));
        
        rentr(ii,2) = rperpP(idx);
        zentr(ii,2) = zperpP(idx);
    
    else % SFD-Minus
        
    end
     
end

if plotit
    plot(rentr(:,1), zentr(:,1), 'or')
    plot(rentr(:,2), zentr(:,2), 'ob')
end

% Compute connection lengths to the divertor target
 
Ldiv = zeros(length(psiSOL),2);
rdiv = zeros(length(psiSOL),2);
zdiv = zeros(length(psiSOL),2);
 
for ii = 1:length(psiSOL)
    
    if psixSL > psixPL % SFD-Plus
        
        % HFS SOL
        
        r0 = rentr(ii,1);
        z0 = zentr(ii,1);
        
        [L, r, z] = calcConnectionLength(r0, z0, psizr257257, ...
            rg257, zg257, bzero, rzero, limdata, -1);
        
        Ldiv(ii,1) = L;
        rdiv(ii,1) = r;
        zdiv(ii,1) = z;
        
        % LFS SOL
        
        r0 = rentr(ii,2);
        z0 = zentr(ii,2);
        
        [L, r, z] = calcConnectionLength(r0, z0, psizr257257, ...
            rg257, zg257, bzero, rzero, limdata, 1);
        
        Ldiv(ii,2) = L;
        rdiv(ii,2) = r;
        zdiv(ii,2) = z;
    
    else % SFD-Minus
        
    end
    
end

if plotit
    plot(rdiv(:,1), zdiv(:,1), 'or')
    plot(rdiv(:,2), zdiv(:,2), 'ob')
end
 
%...............................................................................
% Compute the heat flux profile on the divertor via the diffusion eq

% Heat diffusion coefficient

chi = 0.01; % thermal  diffusivity [Wb^2/s]

% Ion and electron temperatures

Ti = 50 * 11600;  % [K]
Te = 50 * 11600;  % [K]

k = 1.38e-23;  % Boltzmann constant [(m^2*kg)/(s^2*K)]

% Sound speed

mi = 3.344e-27; % Deuterium mass [kg]

cs = sqrt(k*(Ti + Te)/mi);
   
% Compute the time-of-flight for each HFS/LFS fieldline

tauHFS = Ldiv(:,1)/cs;
tauLFS = Ldiv(:,2)/cs;

%.....................................................................
% Compute the heat flux profile at inner SP via the diffusion equation

if psixSL > psixPL % SFD-Plus
    
    idxpsi = 1:length(psiSOL);
    
    minpsi1 = min(psiSOL);
    maxpsi1 = max(psiSOL);
    
    tau = tauHFS;
    
else % SFD-Minus
       
end

psiDivSP1 = linspace(minpsi1, maxpsi1 + 0.005, 100);
  
qdiv_par_SP1 = zeros(length(psiDivSP1),1);
        
for ii = 1:length(idxpsi)-1
    
    idx = idxpsi(ii);
    
    t = tau(idx);
    alpha = sqrt(4*pi*chi*t);
    
    % Define the initial condition
    
    qdiv_OL = q_parallel_IL(idx);
    psiDivFine = linspace(psiSOL(idx), psiSOL(idx+1), 100);
   
   for jj = 1:length(psiDivSP1) % index of points across divertor domain
       
       dpsiDiv = psiDivSP1(jj) - psiDivFine;
       
       int = (qdiv_OL/alpha)*exp(-dpsiDiv.*dpsiDiv/(4*chi*t));
       
       qdiv_par_SP1(jj) = qdiv_par_SP1(jj) + trapz(abs(psiDivFine), int);
       
   end
    
end

rzDiv_SP1 = zeros(length(psiDivSP1),2);
sDiv_SP1  = zeros(length(psiDivSP1),1);

% Determine where the flux surfaces for the SP1 solution intersect the limiter

for ii = 1:length(psiDivSP1)
        
    rz = isoflux_spFinder(psizr257257, psiDivSP1(ii), rg257, zg257, ...
        limdata, limIdxL);
    
    rzDiv_SP1(ii,1) = rz(1,1);
    rzDiv_SP1(ii,2) = rz(1,2);
        
end

for ii = 1:length(psiDivSP1)
    sDiv_SP1(ii) = sLimTot - ...
        calcLimDistance(rzDiv_SP1(ii,1), rzDiv_SP1(ii,2), limdata);  
end

if plotit
    plot(rzDiv_SP1(:,1), rzDiv_SP1(:,2), 'og')
end

% Compute the fieldline incidence angles for SP1 solution

thetaB_SP1 = zeros(length(sDiv_SP1),1);

for ii = 1:length(sDiv_SP1)
    
    r = rzDiv_SP1(ii,1);
    z = rzDiv_SP1(ii,2);
    
    thetaB_SP1(ii) = calcBField_IncidenceAngle(r, z, psizr, rg, zg, ...
        bzero, rzero, limdata);
    
end

qdiv_perp_SP1 = -qdiv_par_SP1.*sin(thetaB_SP1);

%.....................................................................
% Compute the heat flux profile at outer SP via the diffusion equation

if psixSL > psixPL % SFD-Plus
    
    idxps3 = 1:length(psiSOL);
    
    minpsi3 = min(psiSOL);
    maxpsi3 = max(psiSOL);
    
    tau = tauLFS;
    
else % SFD-Minus
    
end

psiDivSP3 = linspace(minpsi3, maxpsi3 + 0.005, 100);
  
qdiv_par_SP3 = zeros(length(psiDivSP3),1);
        
for ii = 1:length(idxpsi)-1
    
    idx = idxpsi(ii);
    
    t = tau(idx);
    alpha = sqrt(4*pi*chi*t);
    
    % Define the initial condition
    
    qdiv_OL = q_parallel_OL(idx);
    psiDivFine = linspace(psiSOL(idx), psiSOL(idx+1), 100);
   
   for jj = 1:length(psiDivSP3) % index of points across divertor domain
       
       dpsiDiv = psiDivSP3(jj) - psiDivFine;
       
       int = (qdiv_OL/alpha)*exp(-dpsiDiv.*dpsiDiv/(4*chi*t));
       
       qdiv_par_SP3(jj) = qdiv_par_SP3(jj) + trapz(psiDivFine, int);
       
   end
    
end

rzDiv_SP3 = zeros(length(psiDivSP3),2);
sDiv_SP3  = zeros(length(psiDivSP3),1);

% Determine where the flux surfaces for the SP3 solution intersect the limiter

for ii = 1:length(psiDivSP3)
        
    rz = isoflux_spFinder(psizr257257, psiDivSP3(ii), rg257, zg257, ...
        limdata, limIdxL);
    
    switch size(rz,1)
        case 2
            rzDiv_SP3(ii,1) = rz(2,1);
            rzDiv_SP3(ii,2) = rz(2,2);
        case 4
            rzDiv_SP3(ii,1) = rz(4,1);
            rzDiv_SP3(ii,2) = rz(4,2);
    end
          
end

for ii = 1:length(psiDivSP3)
    sDiv_SP3(ii) = sLimTot - ...
        calcLimDistance(rzDiv_SP3(ii,1), rzDiv_SP3(ii,2), limdata);  
end

if plotit
    plot(rzDiv_SP3(:,1), rzDiv_SP3(:,2), 'og')
end

% Compute the fieldline incidence angles for SP3 solution

thetaB_SP3 = zeros(length(sDiv_SP3),1);

for ii = 1:length(sDiv_SP3)
    
    r = rzDiv_SP3(ii,1);
    z = rzDiv_SP3(ii,2);
    
    thetaB_SP3(ii) = calcBField_IncidenceAngle(r, z, psizr, rg, zg, ...
        bzero, rzero, limdata);
    
end

qdiv_perp_SP3 = -qdiv_par_SP3.*sin(thetaB_SP3);

if plotit
    
   figure(2)
   hold on
   
   xlabel('s [cm]')
   ylabel('Heat Flux [MW/m^2]')
   title([int2str(shot) ': ' int2str(time) ' ms'])
   
   minQ = min([qdiv_perp_SP1; qdiv_perp_SP3]);
   maxQ = max([qdiv_perp_SP1; qdiv_perp_SP3]);
   
   ylim(1.10*[0 maxQ])
   
   plot(100*[sSPP1 sSPP1], 1.10*[0 maxQ], '-k')
   plot(100*[sSPP2 sSPP2], 1.10*[0 maxQ], '-k')
   
   plot(100*[sSPS1 sSPS1], 1.10*[0 maxQ], '-k')
   plot(100*[sSPS2 sSPS2], 1.10*[0 maxQ], '-k')
   
   plot(100*[s45Deg1 s45Deg1], 1.10*[0 maxQ], '--k')
   plot(100*[s45Deg2 s45Deg2], 1.10*[0 maxQ], '--k')
   
   plot(100*sDiv_SP1, qdiv_perp_SP1, '-r', 'LineWidth', 1)
   plot(100*sDiv_SP3, qdiv_perp_SP3, '-r', 'LineWidth', 1)

end

% Load the IRTV heat flux data
 
qperp = transpose(importdata('qperp_165286_2000.txt'));

qperp = qperp*(100)^2*(1e-6);
 
s = importdata('s_165286_2000.txt');

% Index data for outer SP

idxOuter = find(s > 130);

% Index data for gap

idxGap = find(s < 170);

% Remove the gap

gap1 = s(idxGap(end));
gap2 = s(idxGap(end)+1);

dgap = gap2 - gap1;

s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;

% Plot the data and fit

figure(2)
plot(s, qperp, '-ob', 'LineWidth', 1, 'MarkerSize', 2)

%...............................................................................
% Save the data

HeatFluxSimulation = struct(         ...
    'shot',          shot,           ...
    'time',          time,           ...
    'rg257',         rg257,          ...
    'zg257',         zg257,          ...
    'psizr257257',   psizr257257,    ...
    'rxPL',          rxPL,           ...
    'rxSL',          rxSL,           ...
    'zxPL',          zxPL,           ...
    'zxSL',          zxSL,           ...
    'psixPL',        psixPL,         ...
    'psixSL',        psixSL,         ...
    'rperpP',        rperpP,         ...
    'zperpP',        zperpP,         ...
    'rperpS',        rperpS,         ...
    'zperpS',        zperpS,         ...
    'rentr',         rentr,          ...
    'zentr',         zentr,          ...
    'Ldiv',          Ldiv,           ...
    'rdiv',          rdiv,           ...
    'zdiv',          zdiv,           ...
    'rzDiv_SP1',     rzDiv_SP1,      ...
    'rzDiv_SP3',     rzDiv_SP3,      ...
    'sDiv_SP1',      sDiv_SP1,       ...
    'sDiv_SP3',      sDiv_SP3,       ...
    'qdiv_par_SP1',  qdiv_par_SP1,   ...
    'qdiv_par_SP3',  qdiv_par_SP3,   ...
    'qdiv_perp_SP1', qdiv_perp_SP1', ...
    'qdiv_perp_SP3', qdiv_perp_SP3'  ... 
);

if saveit
    save(['HeatFluxSimulation_' int2str(shot) '_' int2str(time) '.mat'], ...
        'HeatFluxSimulation')
end
