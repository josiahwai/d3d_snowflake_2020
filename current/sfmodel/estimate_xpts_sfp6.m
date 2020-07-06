function xp1 = estimate_xpts_sfp6(eq0, ef, plotit)

c_relax = 0.6;  % relaxation constant on the step sizes

if ~exist('plotit','var'), plotit = 0; end

% analyze equilibrium
struct_to_ws(eq0); clear xlim ylim
snow0 = analyzeSnowflake(eq0);
[rxP0, rxS0] = unpack(snow0.rx);
[zxP0, zxS0] = unpack(snow0.zx);
rSPP = snow0.rSPP([1 end]);  % the 'real' heat flux strike pts
zSPP = snow0.zSPP([1 end]); 
psixP0 = bicubicHermite(rg,zg,psizr,rxP0,zxP0);

% unload eich fit
ef = removeNans(ef);
ef.rsp = ef.rsp([1 end]); % the 'real' heat flux strike pts
ef.zsp = ef.zsp([1 end]); 

% ===================================
% Estimate correction to primary x-pt 
% ===================================

% secondary x-pt should only move to (partially) compensate for strike pt 
% that its closer to. Primary should compensate for the other
dist(1) = norm([rxS0 zxS0] - [rSPP(1) zSPP(1)]);
dist(2) = norm([rxS0 zxS0] - [rSPP(end) zSPP(end)]);
mindist = min(dist); 

% e.g. wtPS = how much Primary x-pt should move to compensate for mismatch
% at Second strike point

% wtPP = 0.9;
% wtSS = 0.4;

if dist(2) < dist(1)
  wtPP = 0.9;
  wtSS = 0.7;
else
  wtPP = 0.9;
  wtSS = 0.3;
end

wtSP = 1-wtPP;
wtPS = 1-wtSS;

% wtPP = 0.9;
% wtSP = 1 - wtPP;
% wtSS = 0.85;
% wtPS = 1 - wtSS;

% Correction due to SP1 
% .....................

% orthogonal projection of strike point error onto flux surfaces
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(1),ef.zsp(1));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
errorP = [ef.rsp(1) - rSPP(1); ef.zsp(1) - zSPP(1)];
resP = u*u'*errorP;



% find (r,z) a distance norm(resP) from (rxP,zxP) and on psibry
th = atan(resP(2)/resP(1));
th_range = linspace(th-pi/5, th+pi/5, 200);
rdum = rxP0 + norm(resP)*cos(th_range);
zdum = zxP0 + norm(resP)*sin(th_range);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
[~,k] = min(abs(psidum - psixP0));

dxpPP = wtPP * [rdum(k) - rxP0; zdum(k) - zxP0];
if dxpPP'*resP < 0, dxpPP = -dxpPP; end  % vector was flipped


% Correction on xp1 due to SP2 mismach
% ....................................
   
% orthogonal projection of strike point error onto flux surfaces
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(2),ef.zsp(2));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
errorS = [ef.rsp(2) - rSPP(2); ef.zsp(2) - zSPP(2)];
resS = u*u'*errorS;

% find (r,z) a distance norm(resP) from (rxP,zxP) and on psibry
th = atan(resS(2)/resS(1));
th_range = linspace(th-pi/5, th+pi/5, 200);
rdum = rxP0 + norm(resS)*cos(th_range);
zdum = zxP0 + norm(resS)*sin(th_range);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
[~,k] = min(abs(psidum - psixP0));

dxpPS = wtPS * [rdum(k) - rxP0; zdum(k) - zxP0];
if dxpPS'*resS < 0, dxpPS = -dxpPS; end  % vector was flipped



% sum contributions together
% ..........................
dxpP = c_relax * (dxpPP + dxpPS); 

rxP1 = rxP0 + dxpP(1);
zxP1 = zxP0 + dxpP(2);

% =====================================
% Estimate correction to secondary x-pt 
% =====================================
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(2),ef.zsp(2));
iso = [dpsidr -dpsidz]' / norm([dpsidr -dpsidz]);
if iso(2) < 0, iso = -iso; end

dxpSS = wtSS * (ef.rsp(2) - rSPP(2)) * iso;

[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(1),ef.zsp(1));
iso = [dpsidr -dpsidz]' / norm([dpsidr -dpsidz]);
if iso(1) < 0, iso = -iso; end

dxpSP = wtSP * (ef.zsp(1) - zSPP(1)) * iso;

dxpS = c_relax * (dxpSS + dxpSP);

rxS1 = rxS0 + dxpS(1);
zxS1 = zxS0 + dxpS(2);


if plotit
  figure(800)
  clf
  plot_eq(eq0,800)
  axis([1 1.5 -1.4 -.9])
  plot(rxS1,zxS1,'bx','markersize',10,'linewidth', 4)
  plot(rxP1,zxP1,'bx','markersize',10,'linewidth', 4)
end

xp1 = [rxP1 rxS1 zxP1 zxS1];

































