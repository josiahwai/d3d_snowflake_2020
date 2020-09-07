% Analyze and return snowflake parameters

function snow = analyzeSnowflake(eq, plotit)

if nargin == 1, plotit = 0; end
if isfield(eq,'gdata'), eq = eq.gdata; end
struct_to_ws(eq);

% find snowflake
[psizr, rg, zg] = regrid(rg, zg, psizr, 257, 257);
[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);

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

[rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 8000, true);
k= find(rlim < 1.36 & rlim > 0.8 & zlim < -1 & zlim > -1.4); % bounding box
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
[~,iS] = findpeaks(-abs(psilim - psixS), 'minpeakheight',-.005);

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

[rSPP, zSPP, sSPP, rSPS, zSPS, sSPS] = removeNans(rSPP, zSPP, sSPP, rSPS, zSPS, sSPS);

if length(sSPS) > 4
  sSPS = sSPS(1:4);
  rSPS = rSPS(1:4);
  zSPS = zSPS(1:4);
end
  

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

snow = struct('rx', [rxP rxS], 'zx',[zxP zxS], 'rSnow', ...
  rSnow, 'zSnow', zSnow, 'drSnow', drSnow, 'dzSnow', dzSnow, 'rSPP', ...
  rSPP, 'zSPP', zSPP, 'sSPP', sSPP, 'rSPS', rSPS, 'zSPS', zSPS, 'sSPS', ...
  sSPS, 'snowType', snowType, 'psixPL',psixP, ....
  'psixSL',psixS, 'snowPlus', snowPlus, 'rg', rg, 'zg', zg, 'psizr', psizr);












