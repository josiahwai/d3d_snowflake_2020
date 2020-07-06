function D = deriv2(x,y,f,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  D = deriv2(x,y,f,n)
%
%  PURPOSE:  Calculate derivatives of f w.r.t. x and y
%
%  INPUTS: x, y = coordinates (need not be equidistant)
%          f = matrix of size(ny,nx) with values on grid x,y
%          n = number of adjacent points to use (default 2)
%
%  OUTPUTS: D = struct matrix with derivatives up to order n
%           D(i,j).f is (j-1):th x, (i-1):th y derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  WRITTEN BY: Anders Welander ON 2017-06-07
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
  n = 2;
end

[ny,nx] = size(f);

if ny ~= numel(y) | nx ~= numel(x)
  error('f must be of size(ny,nx), where ny = numel(y), nx = numel(x)')
end

% Rescale x, y for better numerical conditions
if ny > 1
  dyav = (y(ny)-y(1))/(ny-1);
else
  dyav = 1;
end
if nx > 1
  dxav = (x(nx)-x(1))/(nx-1);
else
  dxav = 1;
end
y = y/dyav;
x = x/dxav;

% Number of points to use for evaluation of derivatives
ndy = min(ny,n+1);
ndx = min(nx,n+1);

% Maximum index difference to different point
di = round((ndy-1)/2);
dj = round((ndx-1)/2);

% Allocate memory for matrices that solve for derivatives
my = zeros(ndy);
mx = zeros(ndx);

% Allocate memory for derivatives
D = repmat(struct('f',f),ndy,ndx);

% Calculate ndy-1 derivatives w.r.t. y
for i = 1:ny
  % Points i1:i2 will be used to evaluate derivatives at i
  i1 = max(1,min(ny-ndy+1,i-di));
  i2 = i1+ndy-1;
  for ie = 1:ndy % index to equation
    DY = y(i1+ie-1)-y(i); % difference in y
    den = 1;
    for id = 1:ndy % index to derivatives
      % f(i1:i2) = my*[0:th to (ndy-1):th derivative]
      my(ie,id) = DY^(id-1)/den;
      den = den*id;
    end
  end
  myi = inv(my);
  h = 1;
  for id = 2:ndy
    h = h/dyav;
    D(id,1).f(i,:) = h*myi(id,:)*f(i1:i2,:);
  end
end

% Calculate ndx-1 derivatives w.r.t. x
for j = 1:nx
  % Points j1:j2 will be used to evaluate derivatives at j
  j1 = max(1,min(nx-ndx+1,j-dj));
  j2 = j1+ndx-1;
  for je = 1:ndx % index to equation
    DX = x(j1+je-1)-x(j); % difference in x
    den = 1;
    for jd = 1:ndx % index to derivatives
      % f(j1:j2) = mx*[0:th to (ndx-1):th derivative]
      mx(je,jd) = DX^(jd-1)/den;
      den = den*jd;
    end
  end
  mxi = inv(mx);
  for id = 1:ndy % Evaluate for all y-derivatives
    h = 1;
    for jd = 2:ndx
      h = h/dxav;
      D(id,jd).f(:,j) = h*D(id,1).f(:,j1:j2)*mxi(jd,:)';
    end
  end
end

