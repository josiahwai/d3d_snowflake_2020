function xp1 = estimate_xpts(eq0,sim0,plotit)
  
if ~exist('plotit','var'), plotit = 0; end
if plotit, figure(800); clf; hold on; end

struct_to_ws(eq0); clear xlim ylim
struct_to_ws(sim0);

[zxP0, zxS0] = unpack(zx);
[rxP0, rxS0] = unpack(rx);
[psixP0, psixS0] = unpack(psix);

% ===================================
% Estimate position of  primary x-pt 
% ===================================

% Correction due to SP1 mismatch
% ...............................

% orthogonal projection of strike point error onto flux surfaces
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(1),z_qmax(1));
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

% Correction due to SP2 mismatch
% ..............................
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(2),z_qmax(2));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
v = [r_qirmax(2) - r_qmax(2); z_qirmax(2) - z_qmax(2)];
res = u*u'*v;

% project a distance |res| orthogonal to previous step
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,rxP0+dxpP1(1), zxP0+dxpP1(2));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
dxpP2 = norm(res)*u;
if sign(dxpP2(1)) ~= sign(res(1)), dxpP2 = -dxpP2; end

dxpP = dxpP1 + dxpP2; 

rxP1 = rxP0 + dxpP(1);
zxP1 = zxP0 + dxpP(2);


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

% find abcd vectors
q = [nesw(end,:); nesw];
abcd = nesw - diff(q) / 2;
abcd = normalize(abcd,2,'norm');


% Determine which vector direction (a,b,c or d)  to move the secondary x-pt
% based on: 1) does the heat flux peak need to move in or out, 2) does x-pt
% need to move to higher or lower psi to match heat flux split

dpsi_abcd = bicubicHermite(rg,zg,psizr,rxS0 + dl*abcd(:,1), zxS0 + dl*abcd(:,2)) - psixS0;
dr_abcd = dl*abcd(:,1);
dr_sp3 = r_qirmax(3) - r_qmax(3);

k = find(sign(dpsi_abcd) == sign(dpsixS) & sign(dr_abcd) == sign(dr_sp3));

abcd_vec = abcd(k,:)';  % direction to move xpt

c_relax = 0.3;
dxp_hfsplit = c_relax * 0.5 * dpsixS / dpsi_abcd(k) * dl * abcd_vec;

s = dxp_hfsplit;
plot(rxS0 + [0 s(1)], zxS0 + [0 s(2)],'g','linewidth',2)

% ========================================
% Secondary x-pt correction: SP3 mismatch 
% ========================================

% SP3 mismatch projection
% ...............................
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,r_qmax(3),z_qmax(3));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
v = [r_qirmax(3) - r_qmax(3); z_qirmax(3) - z_qmax(3)];
r1 = c_relax * u*u'*v;

plot(r_qmax(3) + [0 r1(1)], z_qmax(3) + [0 r1(2)], 'b','linewidth',2)

% option 1

% r = r1;
% rhat = r/norm(r);
% s = dxp_hfsplit;
% dxp_total = s + rhat * (r-s)'*rhat;
%
% plot(rxS0 + [0 r1(1)], zxS0 + [0 r1(2)], 'b','linewidth',2)
% dx = s'*rhat * rhat;
% plot(rxS0 + [0 dx(1)], zxS0 + [0 dx(2)], 'r','linewidth',2)


% option 2

[~,k] = max(abs(nesw*r1));
nesw_vec = nesw(k,:)';
r2 = norm(r1) * sign(nesw_vec'*r1) * nesw_vec;
r = r2;
rhat = r/norm(r);
s = dxp_hfsplit;
dxp_total = s + rhat * (r-s)'*rhat;

plot(rxS0 + [0 r2(1)], zxS0 + [0 r2(2)], 'b','linewidth',2)
dx = s'*rhat * rhat;
plot(rxS0 + [0 dx(1)], zxS0 + [0 dx(2)], 'r','linewidth',2)





% option 3

% find nesw vector most aligned with r1:=residual (should be n/s)
% [~,k] = max(abs(nesw*r1));
% nesw_vec = nesw(k,:)';
% r2 = norm(r1) * sign(nesw_vec'*r1) * nesw_vec;
% r2hat_perp = [r2(2) -r2(1)]' / norm(r2); % unit vector perp to r2
% s = dxp_hfsplit;
% w = s'*(s-r2) / (s'*r2hat_perp) * r2hat_perp;
% dxp_total = r2 + w;


rxS1 = rxS0 + dxp_total(1);
zxS1 = zxS0 + dxp_total(2);


if plotit
  figure(800)
  plot_eq(eq0)
  axis([1 1.5 -1.4 -.9])
  scatter(rxP1,zxP1,'b','filled')
  scatter(rxS1,zxS1,'b','filled')
end

xp1 = [rxP1 rxS1 zxP1 zxS1];















