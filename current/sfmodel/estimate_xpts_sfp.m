function xp1 = estimate_xpts_sfp(eq0, shot, time_ms, plotit)

c_relax = 0.8;  % relaxation constant on the step sizes

if ~exist('plotit','var'), plotit = 0; end

% analyze equilibrium
struct_to_ws(eq0); clear xlim ylim
snow0 = analyzeSnowflake(eq0);
[rxP0, rxS0] = unpack(snow0.rx);
[zxP0, zxS0] = unpack(snow0.zx);
rSPP = snow0.rSPP([1 end]);  % the 'real' heat flux strike pts
zSPP = snow0.zSPP([1 end]); 
psixP0 = bicubicHermite(rg,zg,psizr,rxP0,zxP0);

% analyze heat flux
% .................

% Load heat flux data q(s,t), s=distance along limiter, and t=time
root = '/u/jwai/d3d_snowflake_2020/current/';
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
[~,k] = min(abs(t-time_ms));
qperp = qperp(k,:)';

% obtain parameters from the eich fit to heat flux (strike points etc.)
load('d3d_obj_mks_struct_6565.mat')
ef = eich_fitter(s', qperp, eq0, tok_data_struct);
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
[~,ipkP] = max(dist);
[mindist,ipkS] = min(dist); 

% weights on how much to move the primary x-pt versus secondary x-pt to
% compensate for SP(ipkS) mismatch. Move secondary x-pt more if its closer
wtP = mindist / nansum(dist);           
wtS = 1-wtP;

errorP = [ef.rsp(ipkP) - rSPP(ipkP); ef.zsp(ipkP) - zSPP(ipkP)];
errorS = [ef.rsp(ipkS) - rSPP(ipkS); ef.zsp(ipkS) - zSPP(ipkS)];

correct_secondary_mismatch = 1;
if norm(errorP) > 3*norm(errorS), correct_secondary_mismatch = 0; end  


% Correction due to SP(ipkP) mismatch -- usually SP1
% ..................................................

% orthogonal projection of strike point error onto flux surfaces
[~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(ipkP),ef.zsp(ipkP));
u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
res = u*u'*errorP;



% find (r,z) a distance norm(res) from (rxP,zxP) and on psibry
th = atan(res(2)/res(1));
th_range = linspace(th-pi/5, th+pi/5, 200);
rdum = rxP0 + norm(res)*cos(th_range);
zdum = zxP0 + norm(res)*sin(th_range);
psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
[~,k] = min(abs(psidum - psixP0));

dxpP1 = [rdum(k) - rxP0; zdum(k) - zxP0];
if dxpP1'*res < 0, dxpP1 = -dxpP1; end  % vector was flipped


% Correction on xp1 due to SP(ipkS) mismatch -- usually SP3
% .........................................................
dxpP2 = [0 0];
if correct_secondary_mismatch    
  
  % orthogonal projection of strike point error onto flux surfaces
  [~,dpsidr,dpsidz] = bicubicHermite(rg,zg,psizr,ef.rsp(ipkS),ef.zsp(ipkS));
  u = [dpsidr dpsidz]' / norm([dpsidr dpsidz]);
  resS = u*u'*errorS;
  
  % find (r,z) a distance norm(res) from (rxP,zxP) and on psibry
  th = atan(resS(2)/resS(1));
  th_range = linspace(th-pi/5, th+pi/5, 200);
  rdum = rxP0 + norm(resS)*cos(th_range);
  zdum = zxP0 + norm(resS)*sin(th_range);
  psidum = bicubicHermite(rg,zg,psizr,rdum,zdum);
  [~,k] = min(abs(psidum - psixP0));
  
  dxpP2 = wtP * [rdum(k) - rxP0; zdum(k) - zxP0];
  if dxpP2'*resS < 0, dxpP2 = -dxpP2; end  % vector was flipped
end


% sum contributions together
% ..........................
dxpP = c_relax * (dxpP1 + dxpP2); 

rxP1 = rxP0 + dxpP(1);
zxP1 = zxP0 + dxpP(2);

% =====================================
% Estimate correction to secondary x-pt 
% =====================================
dxpS = [0 0];
if correct_secondary_mismatch  
  %   % orthogonal projection of strike point error onto flux surfaces
  %   dxpS = c_relax * wtS * resS;
  
  
  % isoflux line w
  % ..............
  
  % find the point on divertor leg nearest to secondary x-pt
  crz = contourc(rg,zg,psizr,[psixP0 psixP0]);
  snow0 = analyzeSnowflake(eq0);
  
  rboxmin = min(rxP0, rSPP(ipkS)) - .05;
  rboxmax = max(rxP0, rSPP(ipkS)) + .05;
  zboxmin = min(zxP0, zSPP(ipkS)) - .05;
  zboxmax = max(zxP0, zSPP(ipkS)) + .05;
  
  k = find(crz(1,:) < rboxmax & crz(1,:) > rboxmin & crz(2,:) < zboxmax & ...
    crz(2,:) > zboxmin);
  
  [r,z] = interparc(crz(1,k),crz(2,k), 100, 1);
  iNear = dsearchn([r z], [rxS0 zxS0]);
  
  rNear = r(iNear);  % point on divertor leg nearest secodary x-pt
  zNear = z(iNear);
  
  w = [rxP0 zxP0]' - [rNear zNear]';
  what = w / norm(w);
  
  dxpS = c_relax * wtS * resS'*resS / (what'*resS) * what;
  
  %   scatter(rNear,zNear,'filled')
  %   plot(rNear + [0 what(1)], zNear + [0 what(2)])
  %   plot(rNear + [0 p(1)], zNear + [0 p(2)], 'g','linewidth',2)
  %   plot(rstrike + [0 resS(1)], zstrike + [0 resS(2)], 'g','linewidth',2)
  %   plot(crz(1,k), crz(2,k), 'k','linewidth',2)    
end

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

































