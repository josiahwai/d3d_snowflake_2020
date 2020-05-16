function xp1 = estimate_xpts_sfm(eq0,sim0,plotit)

c_relaxP = 0.8;
c_relaxS = 0.3;  % relaxation constant on the step sizes
thresh = .005;   % don't relax step sizes if under .5 cm

if ~exist('plotit','var'), plotit = 0; end
if plotit, figure(800); clf; hold on; end

struct_to_ws(eq0); clear xlim ylim
struct_to_ws(sim0);

snow0 = analyzeSnowflake(eq0);
[rxP0, rxS0] = unpack(snow0.rx);
[zxP0, zxS0] = unpack(snow0.zx);
rSPP = snow0.rSPP;
zSPP = snow0.zSPP;
psixP0 = bicubicHermite(rg,zg,psizr,rxP0,zxP0);
psixS0 = bicubicHermite(rg,zg,psizr,rxS0,zxS0);


% ===================================
% Estimate position of  primary x-pt 
% ===================================

% Correction due to SP1 mismatch
% ...............................

% orthogonal projection of strike point error onto flux surfaces
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,rSPP(1), zSPP(1));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
v = [r_qirmax(1) - r_qmax(1); z_qirmax(1) - z_qmax(1)];
res = u*u'*v;

% find (r,z) a distance norm(res) from (rxP,zxP) and on psibry
th = atan(res(2)/res(1));
th_range = linspace(th-pi/6, th+pi/6, 200);
rdum = rxP0 + norm(res)*cos(th_range);
zdum = zxP0 + norm(res)*sin(th_range);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
[~,k] = min(abs(psidum - psixP0));

dxpP1 = [rdum(k) - rxP0; zdum(k) - zxP0];
if dxpP1'*res < 0, dxpP1 = -dxpP1; end  % vector was pointed wrong direction


% Correction due to SP2 mismatch
% ..............................
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,rSPP(2),zSPP(2));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
v = [r_qirmax(2) - r_qmax(2); z_qirmax(2) - z_qmax(2)];
res = u*u'*v;

% project a distance |res| orthogonal to previous step
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,rxP0, zxP0+.02);
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
dxpP2 = norm(res)*u;
if sign(dxpP2(1)) ~= sign(res(1)), dxpP2 = -dxpP2; end

dxpP = dxpP1 + dxpP2;

% ==============================================
% Secondary x-pt correction: heat flux splitting 
% ==============================================

% find psixS by determining how heat flux divides between SP2 and SP3
[rmid,k] = max(rbbbs);
zmid = zbbbs(k);
lambda_q = .002; % sol width

drsplit_mid = -lambda_q * log(qirmax(3)/qmax(3) * sum(qmax(2:3))/sum(qirmax(2:3)));

[~,dpsidr] = bicubicHermite(rg,zg,psizr,rmid,zmid);
dpsixS = dpsidr * drsplit_mid;


% find North,East,South,West and A,B,C,D vectors
% ..............................................

% nesw: 4 contour lines at the x-pt, N:=line in the upper right quadrant
% abcd: 4 lines in the +- grad(psi) directions (except that grad(psi) is 0 
%       on the x-pt itself), i.e. halfway between nesw vectors

th = pi/2:-.01:-3/2*pi;
dl = .02;
rdum = rxS0 + dl*cos(th);
zdum = zxS0 + dl*sin(th);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);

k = find((psidum(1:end-1)-psixS0).*(psidum(2:end)-psixS0)<0);  % index nesw

if length(k) == 4  
  nesw = [rdum(k) - rxS0; zdum(k) - zxS0]' / dl;   
else
  warning('trouble finding nesw')
end

fexp = 20;
dxp_split = fexp * drsplit_mid * nesw(2,:)';

if plotit
  plot(rxS0 + [0 dxp_split(1)], zxS0 + [0 dxp_split(2)],'g','linewidth',2)
end

% ========================================
% Secondary x-pt correction: SP3 mismatch 
% ========================================

% SP3 mismatch projection
% ...............................
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(3)-.04,z_qmax(3));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
v = [r_qirmax(3) - r_qmax(3); z_qirmax(3) - z_qmax(3)];
r1 = u*u'*v;  % orthogonal projection of error at strike pt

if plotit
  plot(r_qmax(3) + [0 r1(1)], z_qmax(3) + [0 r1(2)], 'b','linewidth',2)
end

n = nesw(1,:)';
r2 = sign(r1'*n) * norm(r1) * n;
dxpS = dxp_split + r2;

if plotit
  plot(rxS0 + [0 r2(1)], zxS0 + [0 r2(2)], 'b','linewidth',2)
end



% relaxation of the step sizes
if c_relaxS * norm(dxpS) > thresh
  dxpS = c_relaxS*dxpS;
elseif norm(dxpS) > thresh
  dxpS = thresh / norm(dxpS) * dxpS;
end

if c_relaxP * norm(dxpP) > thresh
  dxpP = c_relaxP*dxpP;
elseif norm(dxpP) > thresh
  dxpP = .005 / norm(dxpP) * dxpP;
end

rxP1 = rxP0 + dxpP(1);
zxP1 = zxP0 + dxpP(2);
rxS1 = rxS0 + dxpS(1);
zxS1 = zxS0 + dxpS(2);

if plotit
  figure(800)
  plot_eq(eq0,800)
  axis([1 1.5 -1.4 -.9])
  scatter(rxP1,zxP1,'b','filled')
  scatter(rxS1,zxS1,'b','filled')
end

xp1 = [rxP1 rxS1 zxP1 zxS1];





% x = linspace(1, 1.35, 10);
% y = linspace(-1.45, -1.25, 10);
% [x,y] = meshgrid(x,y);
% f = calc_fluxexp(eqs{1},x,y);
% mean(mean(f))








