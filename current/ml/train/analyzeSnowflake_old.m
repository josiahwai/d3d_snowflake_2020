% Analyze and return snowflake parameters

function snow = analyzeSnowflake(eq, plotit)

if nargin == 1, plotit = 0; end
if isfield(eq,'gdata'), eq = eq.gdata; end
struct_to_ws(eq);

% find snowflake
[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);
[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);

% [rguess,zguess] = isoflux_xpFinder(psizr,1.15,-1.25,rg,zg); 
% [rxP, zxP, rxS, zxS] = snowFinder(psizr, rguess, zguess, 0.1, rg, zg); 
% 
% if rxP < min(rg), rxP = min(rg) + .02; end
% if rxS < min(rg), rxS = min(rg) + .02; end
% if zxP < min(zg), zxP = min(zg) + .02; end
% if zxS < min(zg), zxS = min(zg) + .02; end
% 
% % zoom in on snowflake x-pts
% [rxP, zxP, psixP] = isoflux_xpFinder(psizr, rxP, zxP, rg, zg);
% [rxS, zxS, psixS] = isoflux_xpFinder(psizr, rxS, zxS, rg, zg);
% 
% if abs(psixS - psibry) < abs(psixP - psibry)
%     swap(psixP, psixS);
%     swap(rxP, rxS);
%     swap(zxP, zxS);   
% end

snowPlus=0; snowMinLFS=0; snowMinHFS=0;
if psixS > psixP
  snowPlus = true;
  snowType = 'snowflake plus';
else
  % check the sign of the angle defined by the secondary x-pt, primary x-pt, and
  % plasma axis (with primary x-pt in the middle) to determine HFS/LFS
  % This is more robust than comparing rxP, rxS
  a = [eq.rmaxis eq.zmaxis 0] - [rxP zxP 0];
  b = [rxS zxS 0] - [rxP zxP 0];
  c = cross(a,b);
  
  if c(3) < 0  % CW rotation from a to b
    snowMinLFS = true;
    snowType = 'snowflake minus LFS';
  else
    snowMinHFS = true;
    snowType = 'snowflake minus HFS';
  end
end

% snow parameters
rSnow = (rxP + rxS) / 2;
zSnow = (zxP + zxS) / 2;
drSnow = abs(rxP - rxS);
dzSnow = abs(zxP - zxS);


% Find strike pts and sort according to CCW distance from inner midplane
% this is more robust than isoflux_spFinder.m
% .............................................
tokdir = '/u/jwai/d3d_snowflake_2020/current/inputs/tok_data/';
tokdata = 'd3d_obj_mks_struct_6565.mat';
load([tokdir tokdata]);
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 2000, true);
k= find(rlim < 1.5 & rlim > 0.8 & zlim < -.9 & zlim > -1.4); % bounding box
rlim = rlim(k);
zlim = zlim(k);
slim = sLimTot - calcLimDistance(rlim,zlim,limdata);
psilim = bicubicHermite(rg,zg,psizr,rlim,zlim);

% find and sort primary strike pts
[~,iP] = findpeaks(-abs(psilim - psixP), 'minpeakheight',-.005);

rSPP=[]; zSPP=[]; sSPP=[];
for k = 1:length(iP)
  iNearby = iP(k)-3:iP(k)+3;
  rSPP(k) = interp1(psilim(iNearby), rlim(iNearby), psixP);
  zSPP(k) = interp1(psilim(iNearby), zlim(iNearby), psixP);
  sSPP(k) = interp1(psilim(iNearby), slim(iNearby), psixP);
end
[sSPP,iSort] = sort(sSPP); 
rSPP = rSPP(iSort);
zSPP = zSPP(iSort);

% find and sort secondary strike pts
[~,iS] = findpeaks(-abs(psilim - psixS), 'minpeakheight',-.01);

rSPS=[]; zSPS=[]; sSPS=[];
for k = 1:length(iS)
  iNearby = iS(k)-3:iS(k)+3;
  rSPS(k) = interp1(psilim(iNearby), rlim(iNearby), psixS);
  zSPS(k) = interp1(psilim(iNearby), zlim(iNearby), psixS);
  sSPS(k) = interp1(psilim(iNearby), slim(iNearby), psixS);
end
[sSPS,iSort] = sort(sSPS);
rSPS = rSPS(iSort);
zSPS = zSPS(iSort);


% Find the 'true' heatflux strike pts based on snow geometry
xp2in = inpolygon(rxS, zxS, limdata(2,:), limdata(1,:));
if snowPlus
  sSP_heat = sSPP(1:2);
elseif snowMinLFS
  if xp2in
   sSP_heat = [sSPP(1:2) sSPS([2 4])];
  else
    sSP_heat = sSPP(1:2);
  end
elseif snowMinHFS
  if xp2in
    sSP_heat = [sSPPS([1 3]) sSPP(end-1:end)];
  else
    sSP_heat = sSPP(1:2);
  end
end
sSP_heat = [sSP_heat nan(1,4-length(sSP_heat))];
  

% plot it
% .......
if plotit
  plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
  hold on
  axis equal
  axis([1.0 1.5 -1.4 -0.9])
  contour(rg,zg,psizr,[psixP psixP],'r');
  contour(rg,zg,psizr,[psixS psixS],'b');
  scatter(rSPP,zSPP,'r','filled')
  scatter(rSPS,zSPS,'b','filled')
end

% snow = [rxP rxS zxP zxS rSnow zSnow drSnow dzSnow rSPP zSPP sSPP rSPS zSPS sSPS];
snow = struct('rx', [rxP rxS], 'zx',[zxP zxS], 'rSnow', ...
  rSnow, 'zSnow', zSnow, 'drSnow', drSnow, 'dzSnow', dzSnow, 'rSPP', ...
  rSPP, 'zSPP', zSPP, 'sSPP', sSPP, 'rSPS', rSPS, 'zSPS', zSPS, 'sSPS', ...
  sSPS, 'sSP_heat', sSP_heat, 'snowType', snowType, 'psixPL',psixP, ....
  'psixSL',psixS);












