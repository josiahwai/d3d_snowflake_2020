addpath(genpath('/u/pvail/d3d_snowflake_2019/invMap/functions'))
addpath(genpath('/u/pvail/d3d_snowflake_2019/invMap/gatools'))

fid  = fopen('time.txt', 'r');
time = textscan(fid, '%f');
fclose(fid);
 
time = time{1};
% time = 4800;

shot = 158955;

plotit = 0;
saveit = 1;

shotdir = ['shotdir/' int2str(shot) '/efit03/'];

geqdsk_name = ['g' int2str(shot) '.0' int2str(time)];

%...............................................................................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

limdata = limdata(:,[1:79 81:end]);

% Define strike point segments
 
Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

zShelf = limdata(1,limIdxL(3));

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

efit_dirname = ['/u/pvail/d3d_snowflake_2019/invMap/' shotdir];
    
eq = read_eq(shot, time/1000, efit_dirname);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
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
    contour(rg, zg, psizr, [psixPL psixPL], 'k', 'LineWidth', 2)
    contour(rg, zg, psizr, [psixSL psixSL], 'k', 'LineWidth', 1)
    
    plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)
    plot(rxSL, zxSL, 'xk', 'Markersize', 12, 'LineWidth', 3)
end

% Find the snowflake strike points in the lower divertor

spRZLP = isoflux_spFinder(psizr, psixPL, rg, zg, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr, psixSL, rg, zg, limdata, limIdxL);

% Determine which SPs are the true strike points

spLogic

if plotit
    for ii = 1:size(spRZLP,1)
        plot(spRZLP(ii,1), spRZLP(ii,2), 'ob', 'LineWidth', 2)
    end
    for ii = 1:size(spRZLS,1)
        plot(spRZLS(ii,1), spRZLS(ii,2), 'om', 'LineWidth', 2)
    end
end

% Compute the flux expansion at each strike point

fExpP1 = calcFluxExpansion(spRZLP(1,1), spRZLP(1,2), psizr, rg, zg, psixPL, ...
    bzero, rzero);

fExpP2 = calcFluxExpansion(spRZLP(2,1), spRZLP(2,2), psizr, rg, zg, psixPL, ...
    bzero, rzero);

fExpS1 = calcFluxExpansion(spRZLS(1,1), spRZLS(1,2), psizr, rg, zg, psixPL, ...
    bzero, rzero);

fExpS2 = calcFluxExpansion(spRZLS(2,1), spRZLS(2,2), psizr, rg, zg, psixPL, ...
    bzero, rzero);

% Compute distance to strike points [m]

sSPP1 = sLimTot - calcLimDistance(spRZLP(1,1), spRZLP(1,2), limdata); 
sSPP2 = sLimTot - calcLimDistance(spRZLP(2,1), spRZLP(2,2), limdata);

sSPS1 = sLimTot - calcLimDistance(spRZLS(1,1), spRZLS(1,2), limdata); 
sSPS2 = sLimTot - calcLimDistance(spRZLS(2,1), spRZLS(2,2), limdata);

%...............................................................................
% Modeling of the midplane SOL

% Determine major radius of primary separatrix on the midplane
 
psimid = interp2(rg, zg, psizr, rg, zmaxis);
 
idxP = find(psimid > psixPL, 1, 'last' ); % midplane primary

idxPrimary = [idxP-1 idxP idxP+1 idxP+2];

ppPrimary = spline(rg(idxPrimary), psimid(idxPrimary));

cP = ppPrimary.coefs(2,:);

rootsP = roots([cP(1) cP(2) cP(3) cP(4)-psixPL]) + rg(idxP);

[~,idxminP] = min(abs(rootsP-rg(idxP)));

rmidPrimary = rootsP(idxminP);

% Compute the poloidal and total fields at the midplane

[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, rmidPrimary, zmaxis);

BrMid = -1/(2*pi*rmidPrimary)*dpsidz;
BzMid =  1/(2*pi*rmidPrimary)*dpsidr;

BpMid = sqrt(BrMid*BrMid + BzMid*BzMid);

BTMid = (bzero*rzero)/rmidPrimary;

BTotMid = sqrt(BpMid*BpMid + BTMid*BTMid);

% Radial grid at the plasma edge

lambdaQ = 0.002; % SOL power width [m]
 
rSOLMid = linspace(rmidPrimary + 5e-5, rmidPrimary + 5e-5 + 3*lambdaQ, 100)';

% Compute flux at each grid point

psiSOL = interp2(rg, zg, psizr, rSOLMid, zmaxis);

if plotit
    contour(rg, zg, psizr, [psiSOL(end) psiSOL(end)], '-b')
end

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
    
    rz = isoflux_spFinder(psizr, psiSOL(ii), rg, zg, limdata, limIdxL);
        
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
        
        if rxPL < rxSL % LFS
            
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
        
        else % HFS 
            
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
    
    if 0
        plot(rzLimSOL(ii,1), rzLimSOL(ii,2), 'or')
        plot(rzLimSOL(ii,3), rzLimSOL(ii,4), 'or')
    end
            
end

% Determine the flux tube partitioning

idxInner = find(psiSOL > psixSL);
idxOuter = setdiff(1:length(psiSOL), idxInner);

% Determine the indices of flux tubes which intersect the shelf

idxShelf = find(rzLimSOL(:,4) == zShelf);

% Determine the indices of flux tubes which intersect outer vertical wall

idxVert = find(rzLimSOL(:,3) == limdata(2,74));

% Determine the indices which intersect the bottom limiter

idxBottom = setdiff(1:length(psiSOL), [idxShelf; idxVert]);

%...............................................................................
% Compute the connection lengths for SOL fieldlines

[dpsidr, dpsidz] = rzGrad(psizr, rg, zg);

% Compute curves which are perpendicular to the flux surfaces (primary)

ns = 5e4;

rperpP = zeros(ns+1,1);
zperpP = zeros(ns+1,1);

if psixPL > psixSL
    if rxPL > rxSL % HFS
        rperpP(1) = 1.001*rxPL;
        zperpP(1) = zxPL;
    else           % LFS
        rperpP(1) = rxPL;
        zperpP(1) = 0.999*zxPL;
    end
end

eta0 = 0.01;

for ii = 1:ns
    
    psi_r = interp2(rg, zg, dpsidr, rperpP(ii), zperpP(ii));
    psi_z = interp2(rg, zg, dpsidz, rperpP(ii), zperpP(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpP(ii+1) = rperpP(ii) - eta*psi_r;
    zperpP(ii+1) = zperpP(ii) - eta*psi_z;
    
    psi = interp2(rg, zg, psizr, rperpP(ii+1), zperpP(ii+1));
    if abs((psi - psiSOL(end))) < 1e-5
        break
    end
    
end

rperpP = nonzeros(rperpP);
zperpP = nonzeros(zperpP);

psiperpP = interp2(rg, zg, psizr, rperpP, zperpP);

% Compute curves which are perpendicular to the flux surfaces (secondary)

ns = 5e4;

rperpS = zeros(ns+1,1);
zperpS = zeros(ns+1,1);

rperpS(1) = rxSL;
zperpS(1) = 0.999*zxSL;

eta0 = 0.01;

for ii = 1:ns
    
    psi_r = interp2(rg, zg, dpsidr, rperpS(ii), zperpS(ii));
    psi_z = interp2(rg, zg, dpsidz, rperpS(ii), zperpS(ii));
    
    eta = eta0*exp(-10*sqrt(psi_r*psi_r + psi_z*psi_z));
    
    rperpS(ii+1) = rperpS(ii) - eta*psi_r;
    zperpS(ii+1) = zperpS(ii) - eta*psi_z;
    
    psi = interp2(rg, zg, psizr, rperpS(ii+1), zperpS(ii+1));
    if abs((psi - psiSOL(end))) < 1e-5
        break
    end
    
end

rperpS = nonzeros(rperpS);
zperpS = nonzeros(zperpS);

psiperpS = interp2(rg, zg, psizr, rperpS, zperpS);

if plotit
    plot(rperpP, zperpP, '-b', 'LineWidth', 1)
    plot(rperpS, zperpS, '-b', 'LineWidth', 1)
    
    plot([rxPL rxSL], [zxPL zxSL], '-b', 'LineWidth', 1)
end

% Compute intersection points of flux tubes with the divertor entrance
 
rentr = zeros(length(psiSOL),2);
zentr = zeros(length(psiSOL),2);
 
for ii = 1:length(psiSOL)     
    
    psiSOLi = psiSOL(ii);
    
    if psixSL > psixPL % SFD-Plus
    
    else % SFD-Minus
        
        if rxPL < rxSL % LFS
            
            % HFS SOL
            
            [~,idx] = min(abs(psiSOLi - psiperpP));
            
            rentr(ii,1) = rperpP(idx);
            zentr(ii,1) = zperpP(idx);
            
            % LFS SOL
            
            if ismember(ii, idxBottom)
                
                if ismember(ii, idxInner)
                    
                    rz = isoflux_cpFinder(psizr, psiSOLi, rg, zg, ...
                        [rxPL rxSL zxPL zxSL]);
                    
                    rentr(ii,2) = rz(1);
                    zentr(ii,2) = rz(2);
                    
                else
                    
                    [~,idx] = min(abs(psiSOLi - psiperpS));
                    
                    rentr(ii,2) = rperpS(idx);
                    zentr(ii,2) = zperpS(idx);
                    
                end
                
            else
                
                rentr(ii,2) = NaN;
                zentr(ii,2) = NaN;
                
            end
        
        else % HFS
            
            % HFS SOL
            
            if ismember(ii, idxInner)
                
                rz = isoflux_cpFinder(psizr, psiSOLi, rg, zg, ...
                    [rxPL rxSL zxPL zxSL]);
                
                rentr(ii,1) = rz(1);
                zentr(ii,1) = rz(2);
            
            else
                
                [~,idx] = min(abs(psiSOLi - psiperpS));
                
                rentr(ii,1) = rperpS(idx);
                zentr(ii,1) = zperpS(idx);
            
            end
            
            % LFS SOL
            
            if ismember(ii, idxBottom)
                
                [~,idx] = min(abs(psiSOLi - psiperpP));
                
                rentr(ii,2) = rperpP(idx);
                zentr(ii,2) = zperpP(idx);
            
            else
                
                rentr(ii,2) = NaN;
                zentr(ii,2) = NaN;
            
            end
            
        end
        
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
    
    else % SFD-Minus
            
            % HFS SOL
            
            r0 = rentr(ii,1);
            z0 = zentr(ii,1);
            
            [L, r, z] = calcConnectionLength(r0, z0, psizr, rg, zg, ...
                bzero, rzero, limdata, -1);
            
            Ldiv(ii,1) = L;
            rdiv(ii,1) = r;
            zdiv(ii,1) = z;
            
            % LFS SOL
            
            if ismember(ii, idxBottom)
                
                r0 = rentr(ii,2);
                z0 = zentr(ii,2);
                
                [L, r, z] = calcConnectionLength(r0, z0, psizr, rg, zg, ...
                    bzero, rzero, limdata, 1);
                
                % Fieldline hit vertical wall
                if abs(r - limdata(2,74)) < 1e-4 && ...
                    (z > limdata(1,74) && z < limdata(1,73))
                    Ldiv(ii,2) = NaN;
                    rdiv(ii,2) = NaN;
                    zdiv(ii,2) = NaN;
                else
                    Ldiv(ii,2) = L;
                    rdiv(ii,2) = r;
                    zdiv(ii,2) = z; 
                end
            
            else
                
                Ldiv(ii,2) = NaN;
                rdiv(ii,2) = NaN;
                zdiv(ii,2) = NaN;
            
            end
        
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

%................................................................
% Compute the heat flux profile at SP1 via the diffusion equation

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        
        idxpsi = 1:length(psiSOL);
        
        minpsi1 = min(psiSOL);
        maxpsi1 = max(psiSOL);
        
        tau = tauHFS;
        
    else % HFS
        
        idxpsi = idxOuter;
        
        minpsi1 = min(psiSOL(idxOuter));
        maxpsi1 = max(psiSOL(idxOuter));
        
        tau = tauHFS;
       
    end
    
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
        
    rz = isoflux_spFinder(psizr, psiDivSP1(ii), rg, zg, limdata, limIdxL);
    
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

qdiv_perp_SP1 = qdiv_par_SP1.*sin(thetaB_SP1);

if mean(qdiv_perp_SP1) < 0
    qdiv_perp_SP1 = -qdiv_perp_SP1;
end

%................................................................
% Compute the heat flux profile at SP2 via the diffusion equation

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    idxpsi = idxInner;
    
    minpsi2 = min(psiSOL(idxInner));
    maxpsi2 = max(psiSOL(idxInner));
    
    if rxPL < rxSL % LFS
        tau = tauLFS;
        q_parallel = q_parallel_OL;
    else % HFS
        tau = tauHFS;
        q_parallel = q_parallel_IL;
    end
    
end

psiDivSP2 = linspace(minpsi2 - 0.005, maxpsi2 + 0.005, 100);
  
qdiv_par_SP2 = zeros(length(psiDivSP2),1);
        
for ii = 1:length(idxpsi)-1
    
    idx = idxpsi(ii);
    
    t = tau(idx);
    alpha = sqrt(4*pi*chi*t);
    
    % Define the initial condition
    
    qdiv_OL = q_parallel(idx);
    psiDivFine = linspace(psiSOL(idx), psiSOL(idx+1), 100);
   
   for jj = 1:length(psiDivSP2) % index of points across divertor domain
       
       dpsiDiv = psiDivSP2(jj) - psiDivFine;
       
       int = (qdiv_OL/alpha)*exp(-dpsiDiv.*dpsiDiv/(4*chi*t));
       
       qdiv_par_SP2(jj) = qdiv_par_SP2(jj) + trapz(abs(psiDivFine), int);
       
   end
    
end

rzDiv_SP2 = zeros(length(psiDivSP2),2);
sDiv_SP2  = zeros(length(psiDivSP2),1);

% Determine where the flux surfaces for the SP2 solution intersect the limiter

for ii = 1:length(psiDivSP2)
        
    rz = isoflux_spFinder(psizr, psiDivSP2(ii), rg, zg, limdata, limIdxL);
    
    switch size(rz,1)
        case 2
            rzDiv_SP2(ii,1) = rz(2,1);
            rzDiv_SP2(ii,2) = rz(2,2);
        case 4
            if rxPL < rxSL % LFS
                rzDiv_SP2(ii,1) = rz(2,1);
                rzDiv_SP2(ii,2) = rz(2,2);
            else
                rzDiv_SP2(ii,1) = rz(3,1);
                rzDiv_SP2(ii,2) = rz(3,2);
            end
    end
          
end

for ii = 1:length(psiDivSP2)
    sDiv_SP2(ii) = sLimTot - ...
        calcLimDistance(rzDiv_SP2(ii,1), rzDiv_SP2(ii,2), limdata);  
end

if plotit
    plot(rzDiv_SP2(:,1), rzDiv_SP2(:,2), 'og')
end

% Compute the fieldline incidence angles for SP2 solution

thetaB_SP2 = zeros(length(sDiv_SP2),1);

for ii = 1:length(sDiv_SP2)
    
    r = rzDiv_SP2(ii,1);
    z = rzDiv_SP2(ii,2);
    
    thetaB_SP2(ii) = calcBField_IncidenceAngle(r, z, psizr, rg, zg, ...
        bzero, rzero, limdata);
    
end

qdiv_perp_SP2 = qdiv_par_SP2.*sin(thetaB_SP2);

if mean(qdiv_perp_SP2) < 0
    qdiv_perp_SP2 = -qdiv_perp_SP2;
end

%................................................................
% Compute the heat flux profile at SP3 via the diffusion equation

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        
        idxpsi = idxOuter;
        
        minpsi3 = min(psiSOL(idxOuter)); % WRONGGGGG
        maxpsi3 = max(psiSOL(idxOuter));
        
        tau = tauLFS;
        
    else % HFS
        
        idxpsi = find(~isnan(Ldiv(:,2)));
        
        minpsi3 = min(psiSOL(~isnan(Ldiv(:,2))));
        maxpsi3 = max(psiSOL(~isnan(Ldiv(:,2))));
        
        tau = tauLFS;
       
    end
    
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
       
       qdiv_par_SP3(jj) = qdiv_par_SP3(jj) + trapz(abs(psiDivFine), int);
       
   end
    
end

rzDiv_SP3 = zeros(length(psiDivSP3),2);
sDiv_SP3  = zeros(length(psiDivSP3),1);

% Determine where the flux surfaces for the SP3 solution intersect the limiter

for ii = 1:length(psiDivSP3)
        
    rz = isoflux_spFinder(psizr, psiDivSP3(ii), rg, zg, limdata, limIdxL);
    
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

qdiv_perp_SP3 = qdiv_par_SP3.*sin(thetaB_SP3);

if mean(qdiv_perp_SP3) < 0
    qdiv_perp_SP3 = -qdiv_perp_SP3;
end

if plotit
    
   figure(2)
   hold on
   
   xlabel('s [cm]')
   ylabel('Heat Flux [MW/m^2]')
   title([int2str(shot) ': ' int2str(time) ' ms'])
   
   minQ = min([qdiv_perp_SP1; qdiv_perp_SP2; qdiv_perp_SP3]);
   maxQ = max([qdiv_perp_SP1; qdiv_perp_SP2; qdiv_perp_SP3]);
   
   ylim(1.10*[0 maxQ])
   
   plot(100*[sSPP1 sSPP1], 1.10*[0 maxQ], '-k')
   plot(100*[sSPP2 sSPP2], 1.10*[0 maxQ], '-k')
   
   plot(100*[sSPS1 sSPS1], 1.10*[0 maxQ], '-k')
   plot(100*[sSPS2 sSPS2], 1.10*[0 maxQ], '-k')
   
   plot(100*[s45Deg1 s45Deg1], 1.10*[0 maxQ], '--k')
   plot(100*[s45Deg2 s45Deg2], 1.10*[0 maxQ], '--k')
   
   plot(100*sDiv_SP1, qdiv_perp_SP1, '-ob', 'LineWidth', 1, 'MarkerSize', 2)
   plot(100*sDiv_SP2, qdiv_perp_SP2, '-ob', 'LineWidth', 1, 'MarkerSize', 2)
   plot(100*sDiv_SP3, qdiv_perp_SP3, '-ob', 'LineWidth', 1, 'MarkerSize', 2)

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

%........................
% Eich profile fit to SP1

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        sSP  = sSPP1;
        fExp = fExpP1;
    else % HFS
        sSP  = sSPS1;
        fExp = fExpS1;
    end
    
end

optionsInner.StartPoint = [2*max(qdiv_perp_SP1) 0.01 lambdaQ sSP 0];
   
fitEich_SP1 = fit(sDiv_SP1, qdiv_perp_SP1, ftInner, 'problem', fExp, ...
    optionsInner);
 
% Plot the optimized Eich profile

if plotit
    
    s_Eich = linspace(100,170,1000);
    
    qperp_Eich_SP1 = eichProfileInner(s_Eich, ...
        fitEich_SP1.q0, 100*fitEich_SP1.S, 100*fitEich_SP1.lambdaQ, fExp, ...
        100*fitEich_SP1.s0);
    
    figure(2)
    plot(s_Eich, qperp_Eich_SP1, '-r', 'LineWidth', 1)

end

%........................
% Eich profile fit to SP2

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        sSP  = sSPP2;
        fExp = fExpP2;
    else % HFS
        sSP  = sSPP1;
        fExp = fExpP1;
    end
    
end

optionsInner.StartPoint = [2*max(qdiv_perp_SP2) 0.01 lambdaQ sSP 0];
   
fitEich_SP2 = fit(sDiv_SP2, qdiv_perp_SP2, ftInner, 'problem', fExp, ...
    optionsInner);

% Plot the optimized Eich profile

if plotit
    
    s_Eich = linspace(100,170,1000);
    
    qperp_Eich_SP2 = eichProfileInner(s_Eich, ...
        fitEich_SP2.q0, 100*fitEich_SP2.S, 100*fitEich_SP2.lambdaQ, fExp, ...
        100*fitEich_SP2.s0);
    
    figure(2)
    plot(s_Eich, qperp_Eich_SP2, '-r', 'LineWidth', 1)

end

%........................
% Eich profile fit to SP3

if psixSL > psixPL % SFD-Plus
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        sSP  = sSPS2;
        fExp = fExpS2;
    else % HFS
        sSP  = sSPP2;
        fExp = fExpP2;
    end
    
end

optionsOuter.StartPoint = [2*max(qdiv_perp_SP3) 0.01 lambdaQ sSP 0];
   
fitEich_SP3 = fit(sDiv_SP3, qdiv_perp_SP3, ftOuter, 'problem', fExp, ...
    optionsOuter);

% Plot the optimized Eich profile

if plotit
    
    s_Eich = linspace(100,170,1000);
    
    qperp_Eich_SP3 = eichProfileOuter(s_Eich, ...
        fitEich_SP3.q0, 100*fitEich_SP3.S, 100*fitEich_SP3.lambdaQ, fExp, ...
        100*fitEich_SP3.s0);
    
    figure(2)
    plot(s_Eich, qperp_Eich_SP3, '-r', 'LineWidth', 1)

end
 
%...............................................................................
% Save the data

HeatFluxSimulation = struct(         ...
    'shot',          shot,           ...
    'time',          time,           ...
    'rg257',         rg,             ...
    'zg257',         zg,             ...
    'psizr257257',   psizr,          ...
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
    'rzDiv_SP2',     rzDiv_SP2,      ...
    'rzDiv_SP3',     rzDiv_SP3,      ...
    'sDiv_SP1',      sDiv_SP1,       ...
    'sDiv_SP2',      sDiv_SP2,       ...
    'sDiv_SP3',      sDiv_SP3,       ...
    'qdiv_par_SP1',  qdiv_par_SP1,   ...
    'qdiv_par_SP2',  qdiv_par_SP2,   ...
    'qdiv_par_SP3',  qdiv_par_SP3,   ...
    'qdiv_perp_SP1', qdiv_perp_SP1', ...
    'qdiv_perp_SP2', qdiv_perp_SP2', ...
    'qdiv_perp_SP3', qdiv_perp_SP3', ... 
    'fitEich_SP1',   fitEich_SP1,    ...
    'fitEich_SP2',   fitEich_SP2,    ...
    'fitEich_SP3',   fitEich_SP3,    ...
    'sSPP1',         sSPP1,          ...
    'sSPP2',         sSPP2,          ...
    'sSPS1',         sSPS1,          ...
    'sSPS2',         sSPS2,          ...
    's45Deg1',       s45Deg1,        ...
    's45Deg2',       s45Deg2         ...
);

if saveit
    save(['HeatFluxSimulation_' int2str(shot) '_' int2str(time) '.mat'], ...
        'HeatFluxSimulation')
end
