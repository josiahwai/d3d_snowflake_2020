
% Calculation of thermal diffusivity (chi)

% Tokamak definition and equilibrium 

tok_dir  = '/u/pvail/d3d_snowflake_2019/invMap/';
tok_name = 'd3d_obj_mks_struct_129129.mat';

load([tok_dir tok_name])

limdata = tok_data_struct.limdata;

efit_dirname  = '/u/pvail/d3d_snowflake_2019/invMap/shotdir/165286/efit03/';

eq = read_eq(165286, 3.8, efit_dirname);

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

zmaxis = eq.gdata.zmaxis;

psizr  = eq.gdata.psizr;

[psizr257257, rg257, zg257] = regrid(rg, zg, psizr, 257, 257);

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

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

% Configure plot

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
title([int2str(165286) ': ' int2str(3800) ' ms'])

% Define strike point segments

Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

zShelf = limdata(1,limIdxL(3));

% Analyze the lower Divertor (SFD)

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr257257, rExp, zExp, rhoExp, rg257, zg257);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr257257, rxPL, zxPL, rg257, zg257);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr257257, rxSL, zxSL, rg257, zg257);

contour(rg257, zg257, psizr257257, [psixPL psixPL], 'k', 'LineWidth', 2)
contour(rg257, zg257, psizr257257, [psixSL psixSL], 'k', 'LineWidth', 1)
    
plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)
plot(rxSL, zxSL, 'xk', 'Markersize', 12, 'LineWidth', 3)

% Locate the strike points

spRZLP = isoflux_spFinder(psizr257257, psixPL, rg257, zg257, limdata, limIdxL);
spRZLS = isoflux_spFinder(psizr257257, psixSL, rg257, zg257, limdata, limIdxL);

plot(spRZLP(1,1), spRZLP(1,2), 'om', 'LineWidth', 2)
plot(spRZLP(2,1), spRZLP(2,2), 'om', 'LineWidth', 2)

plot(spRZLS(3,1), spRZLS(3,2), 'ob', 'LineWidth', 2)
plot(spRZLS(4,1), spRZLS(4,2), 'ob', 'LineWidth', 2) 

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

contour(rg257, zg257, psizr257257, [psiSOL(end) psiSOL(end)], '-b')

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

rperpS(1) = rxSL;
zperpS(1) = 0.999*zxSL;

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

plot(rperpP, zperpP, '-b', 'LineWidth', 1)
plot(rperpS, zperpS, '-b', 'LineWidth', 1)
    
plot([rxPL rxSL], [zxPL zxSL], '-b', 'LineWidth', 1)

% Compute intersection points of flux tubes with the divertor entrance
 
rentr = zeros(length(psiSOL),1);
zentr = zeros(length(psiSOL),1);
 
for ii = 1:length(psiSOL)     
    
    psiSOLi = psiSOL(ii);
    
    if psixSL > psixPL % SFD-Plus
    
    else % SFD-Minus
        
        if rxPL < rxSL % LFS
            
            % LFS SOL
            
            if ismember(ii, idxBottom)
                
                if ismember(ii, idxInner)
                    
                    rz = isoflux_cpFinder(psizr257257, psiSOLi, ...
                        rg257, zg257, [rxPL rxSL zxPL zxSL]);
                    
                    rentr(ii,1) = rz(1);
                    zentr(ii,1) = rz(2);
                    
                else
                    
                    [~,idx] = min(abs(psiSOLi - psiperpS));
                    
                    rentr(ii,1) = rperpS(idx);
                    zentr(ii,1) = zperpS(idx);
                    
                end
                
            else
                
                rentr(ii,1) = NaN;
                zentr(ii,1) = NaN;
                
            end
        
        else % HFS
            
        end
        
    end
     
end

plot(rentr(:,1), zentr(:,1), 'ob')

% Compute connection lengths to the divertor target
 
Ldiv = zeros(length(psiSOL),1);
rdiv = zeros(length(psiSOL),1);
zdiv = zeros(length(psiSOL),1);
 
for ii = 1:length(psiSOL)
    
    if psixSL > psixPL % SFD-Plus
    
    else % SFD-Minus
            
            
            % LFS SOL
            
            if ismember(ii, idxBottom)
                
                r0 = rentr(ii,1);
                z0 = zentr(ii,1);
                
                [L, r, z] = calcConnectionLength(r0, z0, psizr257257, ...
                    rg257, zg257, bzero, rzero, limdata, 1);
                
                % Fieldline hit vertical wall
                if abs(r - limdata(2,74)) < 1e-4 && ...
                    (z > limdata(1,74) && z < limdata(1,73))
                    Ldiv(ii,1) = NaN;
                    rdiv(ii,1) = NaN;
                    zdiv(ii,1) = NaN;
                else
                    Ldiv(ii,1) = L;
                    rdiv(ii,1) = r;
                    zdiv(ii,1) = z; 
                end
            
            else
                
                Ldiv(ii,1) = NaN;
                rdiv(ii,1) = NaN;
                zdiv(ii,1) = NaN;
            
            end
        
    end
    
end

plot(rdiv(:,1), zdiv(:,1), 'ob')

% Load the IRTV heat flux data
 
qperp = transpose(importdata('qperp_165286_3800.txt'));
s = importdata('s_165286_3800.txt');

% Fit Eich profile to the data
 
qEich = @(q0, S, lambdaQ, fExp, s0, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
    - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S);
 
ft = fittype(qEich);

options = fitoptions(ft);
options.StartPoint = [75 1 1 5 150];
options.Lower = [0 0 0 0 0];
 
fitEich = fit(s, qperp, ft, options);
 
% Extract the coefficients
 
q0      = fitEich.q0;
S       = fitEich.S;
lambdaQ = fitEich.lambdaQ;
fExp    = fitEich.fExp;
s0      = fitEich.s0;

% Define the optimized Eich profile

s_Eich = linspace(0,300,1000);

qperp_Eich = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (s_Eich - s0)./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (s_Eich - s0)./S);
 
% Plot the data and fit
 
figure(2)
plot(s, qperp, '-ob', 'LineWidth', 1, 'MarkerSize', 2)

hold on
plot(s_Eich, qperp_Eich, '-r', 'LineWidth', 1)
 
xlim([min(s) max(s)])
ylim([0 max(qperp)+5])

plot([s0 s0], [0 max(qperp)+10], '--k')
 
xlabel('Distance s [cm]')
ylabel('Heat Flux [W/cm^2]')

title([int2str(165286) ': ' int2str(3800) ' ms'])
 
% Add fit parameters to plot
 
text(80, 60, ['S = ' num2str(S,'%3.1f') ' cm'], 'FontSize', 12)
 
text(80, 50, ['\lambda_{q} = ' num2str(lambdaQ,'%3.1f') ' cm'], 'FontSize', 12)
 
text(80, 40, ['f_{exp} = ' num2str(fExp,'%3.1f')], 'FontSize', 12)

% Ion and electron temperatures

Ti = 50 * 11600;  % [K]
Te = 50 * 11600;  % [K]

k = 1.38e-23;  % Boltzmann constant [(m^2*kg)/(s^2*K)]

% Sound speed

mi = 3.344e-27; % Deuterium mass [kg]

cs = sqrt(k*(Ti + Te)/mi);

% Compute the time-of-flight for each fieldline (splitting region to divertor)
 
taudiv = Ldiv(idxOuter)./cs;

taudiv = taudiv(~isnan(taudiv));

tau = max(taudiv);

L = 30;
tau = L/cs;

% Calculate the thermal diffusivity in units [m^2/s]

chiTemp = (S/100)^2/tau;

% Compute the flux gradient at SSP
 
[~, dpsidr, dpsidz] = bicubicHermite(rg257, zg257, psizr257257, ...
    spRZLS(4,1), spRZLS(4,2));

gradPsi = sqrt(dpsidr*dpsidr + dpsidz*dpsidz);

% Calculate the thermal diffusivity in units [Wb^2/s]

chi = chiTemp*gradPsi*gradPsi;

fprintf("chi = %3.3f m^2/s\n" , chiTemp)
fprintf("chi = %3.3f Wb^2/s\n", chi)
