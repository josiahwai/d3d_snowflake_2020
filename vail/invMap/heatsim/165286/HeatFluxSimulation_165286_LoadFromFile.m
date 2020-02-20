
shot = 165286;

plotit  = 1;
plotit2 = 0;

shotdir = ['shotdir/' int2str(shot) '/efit03/'];

geqdsk_name = ['g' int2str(shot) '.0' int2str(time)];

% Load data

datafile = ['/u/pvail/d3d_snowflake_2019/invMap/heatsim/165286/' ...
    int2str(time) '/HeatFluxSimulation_165286_' int2str(time) '.mat'];
load(datafile)

%...............................................................................
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

%...............................................................................
% Configure the plots

if plotit
    figure(ii)
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
    
    if rxPL < rxSL % LFS
        
        if size(spRZLP,1) == 3
            spRZLP = spRZLP([1 2],:);
        end
        if size(spRZLP,1) == 4
            spRZLP = spRZLP([1 2],:);
        end
        if size(spRZLS,1) == 4
            spRZLS = spRZLS([3 4],:);
        end
           
    else % HFS        
        
        if size(spRZLP,1) == 3
            spRZLP = spRZLP([2 3],:);
        end
        if size(spRZLP,1) == 4
            spRZLP = spRZLP([3 4],:);
        end
        if size(spRZLP,1) == 5
            spRZLP = spRZLP([4 5],:);
        end
        if (size(spRZLS,1) == 3) || (size(spRZLS,1) == 4)
            spRZLS = spRZLS([1 2],:);
        end
       
    end
   
end

if plotit
    for ii = 1:size(spRZLP,1)
        plot(spRZLP(ii,1), spRZLP(ii,2), 'ob', 'LineWidth', 2)
    end
    for ii = 1:size(spRZLS,1)
        plot(spRZLS(ii,1), spRZLS(ii,2), 'om', 'LineWidth', 2)
    end
end

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

lambdaQ = 0.005; % SOL power width [m]
 
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

Pheat = 10; % [MW]

frad = 0.2; % radiated power fraction
fOBL = 0.7; % fraction of power to outboard divertor

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

% Compute curves which are perpendicular to the flux surfaces (primary)

rperpP = HeatFluxSimulation.rperpP;
zperpP = HeatFluxSimulation.zperpP;

% Compute curves which are perpendicular to the flux surfaces (secondary)

rperpS = HeatFluxSimulation.rperpS;
zperpS = HeatFluxSimulation.zperpS;

if plotit
    plot(rperpP, zperpP, '-b', 'LineWidth', 1)
    plot(rperpS, zperpS, '-b', 'LineWidth', 1)
    
    plot([rxPL rxSL], [zxPL zxSL], '-b', 'LineWidth', 1)
end

% Compute intersection points of flux tubes with the divertor entrance

rentr = HeatFluxSimulation.rentr;
zentr = HeatFluxSimulation.zentr;

if plotit
    plot(rentr(:,1), zentr(:,1), 'or')
    plot(rentr(:,2), zentr(:,2), 'ob')
end

% Compute connection lengths to the divertor target

Ldiv = HeatFluxSimulation.Ldiv;
rdiv = HeatFluxSimulation.rdiv;
zdiv = HeatFluxSimulation.zdiv;

if plotit
    plot(rdiv(:,1), zdiv(:,1), 'or')
    plot(rdiv(:,2), zdiv(:,2), 'ob')
end

%...............................................................................
% Compute the heat flux profile on the divertor via the diffusion eq

%................................................................
% Compute the heat flux profile at SP1 via the diffusion equation

rzDiv_SP1 = HeatFluxSimulation.rzDiv_SP1;

if plotit
    plot(rzDiv_SP1(:,1), rzDiv_SP1(:,2), 'og')
end

sDiv_SP1 = HeatFluxSimulation.sDiv_SP1;

qdiv_perp_SP1 = HeatFluxSimulation.qdiv_perp_SP1';

%................................................................
% Compute the heat flux profile at SP2 via the diffusion equation

rzDiv_SP2 = HeatFluxSimulation.rzDiv_SP2;

if plotit
    plot(rzDiv_SP2(:,1), rzDiv_SP2(:,2), 'og')
end

sDiv_SP2 = HeatFluxSimulation.sDiv_SP2;

qdiv_perp_SP2 = HeatFluxSimulation.qdiv_perp_SP2';

%................................................................
% Compute the heat flux profile at SP3 via the diffusion equation

rzDiv_SP3 = HeatFluxSimulation.rzDiv_SP3;

if plotit
    plot(rzDiv_SP3(:,1), rzDiv_SP3(:,2), 'og')
end

sDiv_SP3 = HeatFluxSimulation.sDiv_SP3;

qdiv_perp_SP3 = HeatFluxSimulation.qdiv_perp_SP3';

sSPP1 = HeatFluxSimulation.sSPP1;
sSPP2 = HeatFluxSimulation.sSPP2;

sSPS1 = HeatFluxSimulation.sSPS1;
sSPS2 = HeatFluxSimulation.sSPS2;

fitEich_SP1 = HeatFluxSimulation.fitEich_SP1;
fitEich_SP2 = HeatFluxSimulation.fitEich_SP2;
fitEich_SP3 = HeatFluxSimulation.fitEich_SP3;

s45Deg1 = HeatFluxSimulation.s45Deg1;
s45Deg2 = HeatFluxSimulation.s45Deg2;

if plotit2
    
   figure(kk)
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
   
   s_Eich = linspace(100,170,1000);
   
   qperp_Eich_SP1 = eichProfileInner(s_Eich, ...
       fitEich_SP1.q0, 100*fitEich_SP1.S, 100*fitEich_SP1.lambdaQ, ...
       fitEich_SP1.fExp, 100*fitEich_SP1.s0);
   
   qperp_Eich_SP2 = eichProfileOuter(s_Eich, ...
       fitEich_SP2.q0, 100*fitEich_SP2.S, 100*fitEich_SP2.lambdaQ, ...
       fitEich_SP2.fExp, 100*fitEich_SP2.s0);
   
   qperp_Eich_SP3 = eichProfileOuter(s_Eich, ...
       fitEich_SP3.q0, 100*fitEich_SP3.S, 100*fitEich_SP3.lambdaQ, ...
       fitEich_SP3.fExp, 100*fitEich_SP3.s0);
   
    plot(s_Eich, qperp_Eich_SP1, '-r', 'LineWidth', 1)
    plot(s_Eich, qperp_Eich_SP2, '-r', 'LineWidth', 1)
    plot(s_Eich, qperp_Eich_SP3, '-r', 'LineWidth', 1)
     
end
