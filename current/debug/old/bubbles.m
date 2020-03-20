% plot the uncertainty bubbles around the x-points


% =======================
% USER/FUNCTION INPUTS
% =======================
ccc
dpsi = .001;
plotit = 1;
shot = 165288;
time_ms = 2900;

root = '/u/jwai/d3d_snowflake_2020/current/';
efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
efit_eq = read_eq(shot, time_ms/1000, efit_dir);
eq = efit_eq.gdata;

% load('init')
% eq = init.gdata;



% ========================================
% ANALYZE THE SNOWFLAKE EQUILIBRIUM (EFIT)
% ========================================

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
ng = 257;
[psizr, rg, zg] = regrid(rg, zg, psizr, ng, ng);
rgg = repmat(rg',ng,1);
zgg = repmat(zg,1,ng);

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

% ==========================
% PLOT THE UNCERTAINTY ZONES 
% ==========================

% 8 pts around grid pt, CW starting from upper left corner
neighbors = [-ng-1 -1 ng-1 ng ng+1 +1 -ng+1 -ng];


%.........................
% zone around primary x-pt
iP = find(psizr<psixPL+dpsi & psizr>psixPL-dpsi); % points within dpsi of x-pt

% eliminate pts w/o at least 4 neighbors
% (i.e., identify only the large zone(s) around the x-pts, not scattered pts)
nMin = 4;
for dum = 1:length(iP)
    nNeighbors = sum(ismember(iP+neighbors,iP),2);
    iP = iP(nNeighbors >= nMin);
    if min(nNeighbors) >= nMin, break; end
end

%...........................
% zone around secondary x-pt
iS = find(psizr<psixSL+dpsi & psizr>psixSL-dpsi); % points within dpsi of x-pt

% eliminate pts w/o at least 4 neighbors
for dum = 1:length(iS)
    nNeighbors = sum(ismember(iS+neighbors,iS),2);
    iS = iS(nNeighbors >= nMin);
    if min(nNeighbors) >= nMin, break; end
end

if psixSL + dpsi > psixPL % one big zone
    i = [iP; iS];
    b = boundary(rgg(i), zgg(i));
    bP = i(b);
    bS = i(b);
else                      % two independent zones
    b = boundary(rgg(iP), zgg(iP));
    bP = iP(b);
    b = boundary(rgg(iS), zgg(iS));
    bS = iS(b);
end
    

if plotit
    plot(rgg(bP),zgg(bP), 'g', 'linewidth', 2)
    plot(rgg(bS),zgg(bS), 'g', 'linewidth', 2)
%     scatter(rgg(iS),zgg(iS),'k','filled')
end

% 
% % ==================================
% % COMPARE TO UNCONSTRAINED (CAKE) EQ
% % ==================================
% load('init')
% eq = init.gdata;
% 
% [psizr, rg, zg] = regrid(eq.rg, eq.zg, eq.psizr, ng, ng);
% 
% % Locate the snowflake in the lower divertor
% rExp   =  1.1500;
% zExp   = -1.2500;
% rhoExp =  0.1000;
% 
% [rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);
% 
% [rxPL, zxPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
% [rxSL, zxSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);
% 
% plot(rxPL, zxPL, 'xr', 'Markersize', 12, 'LineWidth', 3)
% plot(rxSL, zxSL, 'xr', 'Markersize', 12, 'LineWidth', 3)
% 
% 
% % ================================
% % COMPARE TO CONSTRAINED (CAKE) EQ
% % ================================
% load('eq')
% 
% [psizr, rg, zg] = regrid(eq.rg, eq.zg, eq.psizr, ng, ng);
% 
% % Locate the snowflake in the lower divertor
% rExp   =  1.1500;
% zExp   = -1.2500;
% rhoExp =  0.1000;
% 
% [rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);
% 
% [rxPL, zxPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
% [rxSL, zxSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);
% 
% plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)
% plot(rxSL, zxSL, 'xk', 'Markersize', 12, 'LineWidth', 3)

























