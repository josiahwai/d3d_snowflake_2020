% function [rxP, zxP, rxS, zxS, rSnow, zSnow, dR, dZ, rho, theta] = ...
%     snowFinder(psizr, r0, z0, rhoExp, rg, zg, psibry)

% function [rxP, zxP, rxS, zxS, rSnow, zSnow, dR, dZ, rho, theta] = ...
%   snowFinder_oversample(rg, zg, psizr, psibry, r0, z0)

psizr = eqs{end}.psizr;
psibry = eqs{end}.psibry;
rg = eqs{end}.rg;
zg = eqs{end}.zg;
% r0 = 1.25;
% z0 = -1.1;
% 
% dr_grid = 0.4;
% dz_grid = 0.4;
% rmin = max(min(rg)+.05, r0-dr_grid/2);
% zmin = max(min(zg)+.05, z0-dz_grid/2);
% 
% % xp and vp := probe pts in the linearized frame
% xp = linspace(rmin, rmin+dr_grid, 8) - r0;  
% vp = linspace(zmin, zmin+dz_grid, 8) - z0;

[rmin,rmax,zmin,zmax] = unpack([0.95 1.55 -1.45 -0.9]);
r0 = (rmin+rmax)/2;
z0 = (zmin+zmax)/2;
xp = linspace(rmin, rmax, 8) - r0;  
vp = linspace(zmin, zmax, 8) - z0;

c = combvec(xp,vp);
xp = c(1,:); 
vp = c(2,:);
np = length(xp);

% =============================================================
% Solve for the 3rd order flux expansion coefficients, co (see ref)
% co := const, l1, l2, q1, q2, q3, c1, c2, c3, c4
% =============================================================

% equality constraints Aeq*c = beq
Aeq = ...
  [0 -1 0 2*r0 0 2*r0 0 0 0 0; 
  0 0 0 0 0 2*r0 6*r0^2 0 2*r0^2 0;
  0 0 0 0 -2*r0 0 0 2*r0^2 0 6*r0^2];

beq = [0 0 0]';

% probe the plasma at np places 
[dpdx, dpdv] = gradpolyterms(xp,vp);
p = polyterms(xp,vp);

[psip, psip_r, psip_z] = bicubicHermite(rg,zg,psizr,r0+xp, z0+vp);

Ap = [p; dpdx; dpdv];
bp = [psip psip_r psip_z]';


% minimize over c, f(c) = norm(Ap*c-bp) subject to Aeq*c = beq
dum = pinv([Ap'*Ap Aeq'; Aeq zeros(3,3)]) * [Ap'*bp; beq];
co = dum(1:10);   % coefficients
lm = dum(11:end); % lagrange multipliers


% Find the x-pts based on the coefficients
[const, l1, l2, q1, q2, q3, c1, c2, c3, c4] = unpack(co);


xi = (q3*c1 + q2*c4)/(3*(c1*c1 + c4*c4));
zeta  = (q2*c1 - q3*c4)/ (3*(c1*c1 + c4*c4));

P = (l2*c4 - l1*c1)/(3*(c1*c1 + c4*c4)) + xi*xi - zeta*zeta;
Q = (l2*c1 + l1*c4)/(6*(c1*c1 + c4*c4)) + xi*zeta;

xx1 = xi + sqrt(P/2 + sqrt(P^2/4 + Q^2));
xx2 = xi - sqrt(P/2 + sqrt(P^2/4 + Q^2));
vx1 = zeta + sign(Q)*sqrt(-P/2 + sqrt(P^2/4 + Q^2));
vx2 = zeta - sign(Q)*sqrt(-P/2 + sqrt(P^2/4 + Q^2));

[rx1, rx2, zx1, zx2] = unpack([xx1 xx2 vx1 vx2] + [r0 r0 z0 z0]);





xg = linspace(min(xp), max(xp), 30);
vg = linspace(min(vp), max(vp), 30);
psizr2 = zeros(length(vg),length(xg));

for j = 1:length(vg)
  for i = 1:length(xg)
    p = polyterms(xg(i),vg(j));
    psizr2(j,i) = p*co;
  end
end

figure
hold on
contour(xg+r0, vg+z0, psizr2, 100)
scatter([rx1 rx2], [zx1 zx2],'filled')
scatter(r0+xp,z0+vp,'filled')
axis([0.8 1.81 -1.81 -0.8])


% =======================
% FUNCTION: polyterms
% =======================

% Return the polynomial terms in the taylor expansion, p = [x v x^2 x*v ... x*v^2 v^3]
% such that f(x,v) = c*p where c is the coefficients
% c = const, l1, l2, q1, q2, q3, c1, c2, c3, c4
         
% Reference:  D.D. Ryutov, "Local properties of the magnetic field
% in a snowflake divertor" 

function p = polyterms(x,v)
  n = length(x);
  p = zeros(n,10);
  for k = 1:n
    p(k,:) = [1 x(k) v(k) x(k)^2 x(k)*v(k) v(k)^2 x(k)^3 ...
      x(k)^2*v(k) x(k)*v(k)^2 v(k)^3];
  end
end

% =======================
% FUNCTION: gradpolyterms
% =======================

% Return the derivative of the polynomial terms in the taylor expansion
% dpdx and dpdv such that df(x,v)/dx = c*dpdx and df(x,v)/dv = c*dpdv

% see equations 9,10 of D.D. Ryutov, "Local properties of the magnetic field
% in a snowflake divertor" 

function [dpdx, dpdv] = gradpolyterms(x,v)
  n = length(x);
  dpdx = zeros(n,10);
  dpdv = zeros(n,10);

  for k = 1:n
    dpdx(k,:) = [0 1 0 2*x(k) 2*v(k) 0 3*x(k)^2 2*x(k)*v(k) v(k)^2 0];
    dpdv(k,:) = [0 0 1 0 2*x(k) 2*v(k) 0 x(k)^2 2*x(k)*v(k) 3*v(k)^2];
  end
end















