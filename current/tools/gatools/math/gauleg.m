function [x, w] = gauleg(n,x1,x2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: [x, w] = gauleg(n,x1,x2)
%
%  PURPOSE: Compute the Gauss-Legendre quadrature abscissas and weights
%           of n internal points for the integration range [x1, x2]
%
%  INPUTS: n = number of points
%          x1, x2 = lower and upper bound for the integral, defaults -1, +1
%
%  OUTPUTS: x = points to evaluate the function at
%           w = weights for each of the points
%
%  EXAMPLE: [x, w] = gauleg(32, -1, +1);
%           fx = @(x)sin(x).^2;
%           Igl = w*fx(x)
%           Ii = integral(fx, -1, +1) % compare Igl and Ii

%  METHOD:  When n < 369 gauleg uses a code by Erik Olofsson which is based
%           on a numerical recipe that is the fastest up to n equals about 368
%           For higher n the Glaser-Liu-Rokhlin fast algorithm is used since
%           this is faster for large n
%
%  WRITTEN BY: Anders Welander ON 2020-01-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  x1 = -1;
end
if nargin < 3
  x2 = +1;
end

if n < 369
  x = ones(n,1);
  w = ones(1,n);
  m = ceil(n/2);
  xm = (x2 + x1)/2;
  xl = (x2 - x1)/2;

  for i = 1:m
    z = cos(pi*(i-.25)/(n+.5));
    z1 = Inf;
    while abs(z-z1) > eps
      p1 = 1;
      p2 = 0;
      for j = 1:n
	p3 = p2;
	p2 = p1;
	p1 = ((2*j-1)*z*p2 - (j-1)*p3)/j;
      end
      pp = n*(z*p1-p2)/(z*z-1);
      z1 = z;
      z = z1 - p1/pp;
    end
    x(i) = xm - xl*z;
    x(1+n-i) = xm + xl*z;
    w(i) = 2*xl/((1-z*z)*pp*pp);
    w(1+n-i) = w(i);
  end

else

  [x ders] = alg0_Leg(n);               % Nodes and P_n'(x)
  w = 2./((1-x.^2).*ders.^2)';          % Quadrature weights
  v = 1./ders; v = v./max(abs(v));      % Barycentric weights
  if ~mod(n,2)
    ii = (floor(n/2)+1):n;
    v(ii) = -v(ii);
  end

  % Normalise so that sum(w) = 2
  w = (2/sum(w))*w;           

  % Rescale to arbitrary finite interval
  if x1 ~= -1 | x2 ~= +1
   dab = x2-x1;
   x = (x+1)/2*dab + x1;
   w = dab*w/2;
  end

end

% -------------------- Routines for FAST algorithm ------------------------

function [roots ders] = alg0_Leg(n) % Driver for 'Fast'.

% Compute coefficients of P_m(0), m = 0,..,N via recurrence relation.
Pm2 = 0;
Pm1 = 1; 
for k = 0:n-1
  P = -k*Pm2/(k+1);
  Pm2 = Pm1;
  Pm1 = P;
end

% Get the first roots and derivative values to initialise.
roots = zeros(n,1);
ders = zeros(n,1); % Allocate storage
if mod(n,2) % n is odd
  roots((n-1)/2) = 0; % Zero is a root
  ders((n+1)/2) = n*Pm2; % P'(0)    
else % n is even
  [roots(n/2+1) ders(n/2+1)] = alg2_Leg(P,n); % Find first root
end       

[roots ders] = alg1_Leg(roots,ders); % Other roots and derivatives

% -------------------------------------------------------------------------

function [roots ders] = alg1_Leg(roots,ders)  % Main algorithm for 'Fast'
n = length(roots);
if mod(n,2)
  N = (n-1)/2;
  s = 1;
else
  N = n/2;
  s = 0;
end   

% Approximate roots via asymptotic formula.
k = (n-2+s)/2:-1:1;
theta = pi*(4*k-1)/(4*n+2);
roots(((n+4-s)/2):end) = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);
x = roots(N+1);

% Number of terms in Taylor expansion.
m = 30;

% Storage
hh1 = ones(m+1,1); zz = zeros(m,1); u = zeros(1,m+1); up = zeros(1,m+1);

% Loop over all the roots we want to find (using symmetry).
for j = N+1:n-1
  % Distance to initial approx for next root (from asymptotic foruma).
  h = roots(j+1) - x;

  % Recurrence Taylor coefficients (scaled & incl factorial terms).
  M = 1/h;                           % Scaling
  c1 = 2*x/M; c2 = 1./(1-x^2);       % Some constants
  % Note, terms are flipped for more accuracy in inner product calculation.
  u([m+1 m]) = [0 ders(j)/M];  up(m+1) = u(m);
  for k = 0:m-2
    up(m-k) = (c1*(k+1)*u(m-k)+(k-n*(n+1)/(k+1))*u(m-k+1)/M^2)*c2;
    u(m-(k+1)) = up(m-k)/(k+2);
  end
  up(1) = 0;  

  % Newton iteration
  hh = hh1; step = inf;  l = 0; 
  while (abs(step) > eps) && (l < 10)
    l = l + 1;
    step = (u*hh)/(up*hh)/M;
    h = h - step;        
    Mhzz = (M*h)+zz;
    hh = [1;cumprod(Mhzz)];     % Powers of h (This is the fastest way!)
    hh = hh(end:-1:1);          % Flip for more accuracy in inner product 
  end

  % Update
  x = x + h;
  roots(j+1) = x;
  ders(j+1) = M*(up*hh);  

end

% Nodes are symmetric.
roots(1:N+s) = -roots(n:-1:N+1);
ders(1:N+s) = ders(n:-1:N+1);

% -------------------------------------------------------------------------

function [x1 d1] = alg2_Leg(Pn0,n) % Find the first root (note P_n'(0)==0)
% Approximate first root via asymptotic formula
k = ceil(n/2);
theta = pi*(4*k-1)/(4*n+2);
x1 = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

m = 30; % Number of terms in Taylor expansion.

% Recurrence Taylor coefficients (scaled & incl factorial terms).
M = 1/x1; % Scaling
zz = zeros(m,1); u = [Pn0 zeros(1,m)]; up = zeros(1,m+1); % Allocate storage
for k = 0:2:m-2
  up(k+2) = (k-n*(n+1)/(k+1))*u(k+1)/M^2;
  u(k+3) = up(k+2)/(k+2);
end
% Flip for more accuracy in inner product calculation.
u = u(m+1:-1:1);
up = up(m+1:-1:1);

% Newton iteration
x1k = ones(m+1,1);
step = inf;
l = 0;
while abs(step) > eps & l < 10
  l = l + 1;
  step = (u*x1k)/(up*x1k)/M;
  x1 = x1 - step;
  x1k = [1;cumprod(M*x1+zz)]; % Powers of h (This is the fastest way!)
  x1k = x1k(end:-1:1);        % Flip for more accuracy in inner product
end

% Get the derivative at this root, i.e. P'(x1).
d1 = M*(up*x1k);


