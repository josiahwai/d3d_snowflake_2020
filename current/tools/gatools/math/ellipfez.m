function [F,E,Z] = ellipfez(u,m)
%   [F,E,Z] = ELLIPFEZ(U,M) returns the values of the incomplete
%   elliptic integrals of the first and second kind and Jacobi's
%   Zeta function for arguments U and M. The arrays U and M must
%   be the same size (or either can be scalar).
%   As currently implemented, M is limited to -inf <= M <= 1. 
%
%   Some definitions of the elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
% 
%   See also ELLIPKE, ELLIPJ

%   ELLIPFEZ uses the method of the Arithmetic-Geometric Mean 
%   and Descending Landen Transformation described in [1] Ch. 17.6,
%   to determine the value of the Incomplete Elliptic Integrals 
%   of the First, Second Kind and Jacobi's Zeta Function [1], [2].
%
%       F(phi,m) = integral(1/sqrt(1-m*sin(t)^2), t=0..phi);
%       E(phi,m) = integral(sqrt(1-m*sin(t)^2), t=0..phi);
%       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989
%
%  VERSION @(#)ellipfez.m	1.1 02/14/12
%
%  WRITTEN BY:  Anders Welander  ON	5/21/11
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin<2, error('Not enough input arguments.'); end
  if any(m > 1), error('M must be in the range -inf < M <= 1.'); end
  if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real.')
  end

  if length(m)==1, m = m+0*u; end
  if length(u)==1, u = u+0*m; end
  
  % Evaluate integral(1/sqrt(a^2*cos(t)^2+b^2*sin(t)^2), t=0..phi)
  phi = u;
  a = 1+0*m;
  b = sqrt(1-m);
  c = sqrt(m);
  e = 1;
  C = 0*a;
  Z = 0*a; % Will hold the Jacobi Zeta function
  
  % Substitute variables and reevaluate limit iteratively
  for j = 1:10
    phi = phi + atan(b./a.*tan(phi))+pi*round(phi/pi);
    C = C + e*c.^2;
    c = (a-b)/2;
    Z = Z+c.*sin(phi);
    tmp = (a+b)/2;
    b = sqrt(a.*b);
    a = tmp;  
    e = e*2;
  end
  % At this point a=b=agm(a,b) so the integrand is 1/a

  F = phi./a/e; % Incomplete elliptic integral of the first kind
  E = (Z+(1-C/2).*F); % Incomplete elliptic integral of the second kind
  
  F(m<0 & isinf(m)) = 0;
  F(m==0) = u(m==0);
  E(m==0) = u(m==0);
  Z(m==0) = 0;

  m1 = find(m == 1);
  um1 = abs(u(m1)); 
  if ~isempty(m1), 
      N = floor( (um1+pi/2)/pi );  
      M = find(um1 < pi/2);              

      F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));   
      F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));

      E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1)); 

      Z(m1) = (-1).^N .* sin(u(m1));                      
  end
