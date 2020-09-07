function xpnew = estimate_xpts_sfp7(eq0, ef, plotit)

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
psixS0 = bicubicHermite(rg,zg,psizr,rxS0,zxS0);

% unload eich fit
ef = removeNans(ef);
ef.rsp = ef.rsp([1 end]); % the 'real' heat flux strike pts
ef.zsp = ef.zsp([1 end]); 

% ===========================
% Assign weights to movementes
% ===========================

% secondary x-pt should only move to (partially) compensate for strike pt 
% that its closer to. Primary should compensate for the other
dist(1) = norm([rxS0 zxS0] - [rSPP(1) zSPP(1)]);
dist(2) = norm([rxS0 zxS0] - [rSPP(end) zSPP(end)]);
mindist = min(dist); 

% weights on movement, e.g., wtP1 := how much Primary x-pt should move to 
% compensate for mismatch at strike point 1

% secondary x-point is closer to strike point 2 (normal)
if dist(2) <= dist(1)
  wtP1 = 0.9;
  wtS1 = 0.1;
  wtP2 = 0.1;
  wtS2 = 0.9;
  
%   wtP1 = 0.5;
%   wtS1 = 0.5;
%   wtP2 = 0.5;
%   wtS2 = 0.5;
end

% secondary x-point is closer to strike point 1
if dist(2) > dist(1)
  wtP1 = 0.8;
  wtS1 = 0.2;
  wtP2 = 0.8;
  wtS2 = 0.2;
end

% ==========================
% Find vectors and geometry
% ==========================

% strike pt unit vectors
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr, rSPP(1), zSPP(1));
gradpsihat1 = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
psi_tilda1 = [-gradpsihat1(2) gradpsihat1(1)]';

[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr, rSPP(2), zSPP(2));
gradpsihat2 = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
psi_tilda2 = [gradpsihat2(2) -gradpsihat2(1)]';

% N and E unit vectors at x-pts
th = pi/4:.01:3/4*pi;
dl = .02;

rdum = rxP0 + dl*cos(th);
zdum = zxP0 + dl*sin(th);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
k = find((psidum(1:end-1)-psixP0).*(psidum(2:end)-psixP0)<0);  
N1 = [rdum(k)-rxP0; zdum(k) - zxP0] / dl;
E1 = [N1(2) -N1(1)]';

rdum = rxS0 + dl*cos(th);
zdum = zxS0 + dl*sin(th);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
k = find((psidum(1:end-1)-psixS0).*(psidum(2:end)-psixS0)<0);  
N2 = [rdum(k)-rxS0; zdum(k) - zxS0] / dl;
E2 = [N2(2) -N2(1)]';


% error vectors for strike points
dsp1 = [ef.rsp(1) - rSPP(1); ef.zsp(1) - zSPP(1)];
dsp2 = [ef.rsp(2) - rSPP(2); ef.zsp(2) - zSPP(2)];

err1 =  gradpsihat1 * gradpsihat1' * dsp1;
err2 =  gradpsihat2 * gradpsihat2' * dsp2;


% ===================
% x-point corrections
% ===================
dxp1 = wtP1 * sign(err1'*N1) * norm(err1) * N1 + ...
       wtP2 * sign(err2'*E1) * norm(err2) * E1;


% xp2 is closer to sp2 (normal)
dxp2 = wtS2 * norm(err2) * psi_tilda2 * sign(err2(1)) + ...
       wtS1 * norm(err1) * psi_tilda1' * gradpsihat2 * gradpsihat2;


xp0 = [rxP0 rxS0 zxP0 zxS0];
dxp = c_relax * [dxp1(1) dxp2(1) dxp1(2) dxp2(2)];
xpnew = xp0 + dxp;

if plotit
  figure(800)
  clf
  plot_eq(eq0,800)
  axis([1 1.5 -1.4 -.9])
  scatter( xpnew(1:2), xpnew(3:4), 100, 'bx')
end

































