function [Acell, RAcell, ARcell] = polycellint(rg,zg,rp,zp)
%
%  USAGE:   [Acell, RAcell, ARcell] = polycellint(rg,zg,rp,zp)
%
%  PURPOSE: Calculate how much area of each cell in grid rg, zg
%           is inside the polygon rp, zp
%
%  INPUTS: rg, zg, coordinates for centers of rectangles
%          rp, zp, coordinates defining a polygon
%
%  OUTPUTS: Acell, area of cells that are inside the polygon
%	    RAcell, surface integral of R over inside area
%           ARcell, surface integral of 1/R over inside area
%           All output sizes are [length(zg),length(rg)]
%
	
%  VERSION @(#)polycellint.m	1.2 02/22/15
%
%  WRITTEN BY:  Anders Welander  ON	7/5/14
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid quantities
nr = length(rg);
nz = length(zg);
ngg = nz*nr;
dr = (rg(nr)-rg(1))/(nr-1);
dz = (zg(nz)-zg(1))/(nz-1);
rgg = ones(nz,1)*rg(:)';
zgg = zg(:)*ones(1,nr);
Ag = dr*dz;
RA = rgg*Ag; % Surface integral of R over cell
AR = (log(rgg+dr/2)-log(rgg-dr/2))*dz; % Surface integral of 1/R over cell

% Polygon quantities
if rp(1) ~= rp(end) && zp(1) ~= zp(end)
  rp(end+1) = rp(1);
  zp(end+1) = zp(1);
end
np = length(rp);
ap = angle(rp-mean(rg)+1i*zp-mean(zg));
dap = diff(ap);
if sum(dap<0) > sum(dap>0)
  rp = rp(np:-1:1);
  zp = zp(np:-1:1);
end
drp = diff(rp);
dzp = diff(zp);

% grid number units, starting at zero
irp = (rp(:)-rg(1))/dr;
izp = (zp(:)-zg(1))/dz;

ir1 = round(min([irp(1:end-1)'; irp(2:end)']));
ir2 = round(max([irp(1:end-1)'; irp(2:end)']))-1;
iz1 = round(min([izp(1:end-1)'; izp(2:end)']));
iz2 = round(max([izp(1:end-1)'; izp(2:end)']))-1;
nedge = 0;
xedge = [];
fedge = [];
gedge = [];
iedge = [];
redge = [];
zedge = [];
for j = 1:np-1
  i1 = nedge+1;
  for k = ir1(j) : ir2(j)
    nedge = nedge+1;
    redge(nedge) = (k+0.5)*dr+rg(1);
    iedge(nedge) = j;
    xedge(nedge) = (redge(nedge)-rp(j))/drp(j);
    zedge(nedge) = zp(j)+xedge(nedge)*dzp(j);
    if drp(j) > 0
      fedge(nedge) = nz;
      gedge(nedge) = (k+1)*nz+round((zedge(nedge)-zg(1))/dz)+1;
    else
      fedge(nedge) = -nz;
      gedge(nedge) = k*nz+round((zedge(nedge)-zg(1))/dz)+1;
    end
  end
  for k = iz1(j) : iz2(j)
    nedge = nedge+1;
    zedge(nedge) = (k+0.5)*dz+zg(1);
    iedge(nedge) = j;
    xedge(nedge) = (zedge(nedge)-zp(j))/dzp(j);
    redge(nedge) = rp(j)+xedge(nedge)*drp(j);
    if dzp(j) > 0
      fedge(nedge) = 1;
      gedge(nedge) = round((redge(nedge)-rg(1))/dr)*nz+k+2;
    else
      fedge(nedge) = -1;
      gedge(nedge) = round((redge(nedge)-rg(1))/dr)*nz+k+1;
    end
  end
  i2 = nedge;
  n = i2-i1+1;
  [xedge(i1:i2), kk(1:n)] = sort(xedge(i1:i2));
  fedge(i1:i2) = fedge(i1-1+kk(1:n));
  gedge(i1:i2) = gedge(i1-1+kk(1:n));
  iedge(i1:i2) = iedge(i1-1+kk(1:n));
  redge(i1:i2) = redge(i1-1+kk(1:n));
  zedge(i1:i2) = zedge(i1-1+kk(1:n));
end
  
redge(nedge+1) =   redge(1);
zedge(nedge+1) =   zedge(1);
iedge(nedge+1) =   iedge(1);
fedge(nedge+1) =   fedge(1);
gedge(nedge+1) =   gedge(1);
xedge(nedge+1) =   xedge(1);

ncell = zeros(nz,nr); % Number of entry and exit points
fcell = zeros(ngg,6); % How grid index (ig) changes when exiting cell, going ccw
rcell = zeros(ngg,6); % The r of entry and exit points
zcell = zeros(ngg,6); % The z of entry and exit points

% Covered area of grid cell
Acell = zeros(nz,nr);

% Surface integral of R over covered area within grid cell
RAcell = zeros(nz,nr);

% Surface integral of 1/R over covered area within grid cell
ARcell = zeros(nz,nr);

for j = 1:nedge
  ig = gedge(j); % Entering grid cell ig:
  if ncell(ig) > 0 % Calculate areas from previous exit (a) to new entry (b)
    ri = rgg(ig)-dr/2;
    ro = rgg(ig)+dr/2;
    zd = zgg(ig)-dz/2;
    zu = zgg(ig)+dz/2;
    za = zcell(ig,ncell(ig));
    zb = zedge(j);
    fa = fcell(ig,ncell(ig));
    fb = fedge(j);
    if fa == -nz % contour went out through inner edge of the cell
      if fb == nz % contour came in through inner edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ri;
        RAcell(ig) = RAcell(ig) + (zb-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ri);
      else % Lower inner corner must be inside polygon
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1-nz) == 0
	  Acell(ig-1-nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zd-za)*ri;
        RAcell(ig) = RAcell(ig) + (zd-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zd-za)*log(ri);
	if fb == -nz % contour came in from outer edge
	   Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	  RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	elseif fb == -1 % contour came in from top
	   Acell(ig) =  Acell(ig) + dz*ro;
	  RAcell(ig) = RAcell(ig) + dz*ro^2/2;
	  ARcell(ig) = ARcell(ig) + dz*log(ro);
	end
      end
    elseif fa == nz % contour went out through outer edge of the cell
      if fb == -nz % contour came in through outer edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ro;
        RAcell(ig) = RAcell(ig) + (zb-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ro);
      else % Upper outer corner must be inside polygon
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1+nz) == 0
	  Acell(ig+1+nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zu-za)*ro;
        RAcell(ig) = RAcell(ig) + (zu-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zu-za)*log(ro);
	if fb == nz % contour came in from inner edge
	   Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	  RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	elseif fb ==  1
	   Acell(ig) =  Acell(ig) - dz*ri;
	  RAcell(ig) = RAcell(ig) - dz*ri^2/2;
	  ARcell(ig) = ARcell(ig) - dz*log(ri);
	end
      end
    elseif fa == -1 % contour went out through bottom of the cell
      if fb ~= 1 % contour didn't come in through bottom of the cell
	% Lower outer corner must be inside polygon
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1+nz) == 0
	  Acell(ig-1+nz) = Ag;
	end
	if fb == -nz % contour came in through outer edge of the cell
           Acell(ig) =  Acell(ig) + (zb-zd)*ro;
          RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	else
           Acell(ig) =  Acell(ig) + dz*ro;
          RAcell(ig) = RAcell(ig) + dz*ro^2/2;
          ARcell(ig) = ARcell(ig) + dz*log(ro);
	  if fb == nz % contour came in from inner edge
             Acell(ig) =  Acell(ig) + (zb-zu)*ri;
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	  end
	end
      end
    elseif fa == 1 % contour went out through top of the cell
      if fb ~= -1 % contour didn't come in through top of the cell
	% Upper inner corner must be inside polygon
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1-nz) == 0
	  Acell(ig+1-nz) = Ag;
	end
	if fb == nz % contour came in from inner edge
           Acell(ig) =  Acell(ig) + (zb-zu)*ri;
          RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	else
           Acell(ig) =  Acell(ig) - dz*ri;
          RAcell(ig) = RAcell(ig) - dz*ri^2/2;
          ARcell(ig) = ARcell(ig) - dz*log(ri);
	  if fb == -nz % contour came in through outer edge of the cell
             Acell(ig) =  Acell(ig) + (zb-zd)*ro;
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	  end
	end
      end
    end    
  end
  ncell(ig) = ncell(ig)+1;
  rcell(ig,ncell(ig)) = redge(j);
  zcell(ig,ncell(ig)) = zedge(j);
  fcell(ig,ncell(ig)) = fedge(j);
  xa = xedge(j);
  k = iedge(j); % k is index to bbbs point
  keep_going = true;
  while keep_going % scan iedge(j):iedge(j+1)
    if k == iedge(j+1)
      xb = xedge(j+1);
    else
      xb = 1;
    end
    % The coefficients of r(x)*dz(x):
    xx(1) = rp(k)*dzp(k);
    xx(2) = drp(k)*dzp(k);
    ya = 1;
    yb = 1;
    dA = 0;
    for i = 1:2
      ya = ya*xa;
      yb = yb*xb;
      dA = dA + xx(i)*(yb-ya)/i;
    end
    Acell(ig) = Acell(ig) + dA;
    % The coefficients of r(x)*r(x)*dz(x):
    x1(3) = drp(k)*xx(2);
    x1(2) = drp(k)*xx(1)+rp(k)*xx(2);
    x1(1) = rp(k)*xx(1);
    ya = 1;
    yb = 1;
    dRA = 0;
    for i = 1:3
      ya = ya*xa;
      yb = yb*xb;
      dRA = dRA + x1(i)*(yb-ya)/i/2;
    end
    RAcell(ig) = RAcell(ig) + dRA;
    if drp(k) ~= 0
      ra = rp(k)+drp(k)*xa;
      rb = rp(k)+drp(k)*xb;
      dAR = dzp(k)*(rb*log(rb)-ra*log(ra)-drp(k)*(xb-xa))/drp(k);
    else
      dAR = dzp(k)*log(rp(k))*(xb-xa);
    end
    ARcell(ig) = ARcell(ig) + dAR;
    xa = 0; % If the k loop continues with next bbbs point then xa will be 0
    if k == iedge(j+1)
      keep_going = false;
    else
      k = k+1;
      if k == np
	k = 1;
      end
    end
  end
  % Exiting cell ig:
  ncell(ig) = ncell(ig)+1;
  rcell(ig,ncell(ig)) = redge(j+1);
  zcell(ig,ncell(ig)) = zedge(j+1);
  fcell(ig,ncell(ig)) = fedge(j+1);
end


% Walk ccw along the edges of the cells from the exit point (a) to the entry (b)
for ig = 1:ngg
  if ncell(ig) > 0
    ri = rgg(ig)-dr/2;
    ro = rgg(ig)+dr/2;
    zd = zgg(ig)-dz/2;
    zu = zgg(ig)+dz/2;
    za = zcell(ig,ncell(ig));
    zb = zcell(ig,1);
    fa = fcell(ig,ncell(ig));
    fb = fcell(ig,1);
    if fa == -nz % contour went out through inner edge of the cell
      if fb == nz % contour came in through inner edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ri;
        RAcell(ig) = RAcell(ig) + (zb-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ri);
	if zb > za % In this case we just calculated what was nicked out
           Acell(ig) =  Acell(ig) + Ag;
          RAcell(ig) = RAcell(ig) + RA(ig);
          ARcell(ig) = ARcell(ig) + AR(ig);	  
	end
      else % Lower inner corner must be inside polygon
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1-nz) == 0
	  Acell(ig-1-nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zd-za)*ri;
        RAcell(ig) = RAcell(ig) + (zd-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zd-za)*log(ri);
	if fb == -nz % contour came in from outer edge
	   Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	  RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	elseif fb == -1 % contour came in from top
	   Acell(ig) =  Acell(ig) + dz*ro;
	  RAcell(ig) = RAcell(ig) + dz*ro^2/2;
	  ARcell(ig) = ARcell(ig) + dz*log(ro);
	end
      end
    elseif fa == nz % contour went out through outer edge of the cell
      if fb == -nz % contour came in through outer edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ro;
        RAcell(ig) = RAcell(ig) + (zb-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ro);
      else % Upper outer corner must be inside polygon
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1+nz) == 0
	  Acell(ig+1+nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zu-za)*ro;
        RAcell(ig) = RAcell(ig) + (zu-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zu-za)*log(ro);
	if fb == nz % contour came in from inner edge
	   Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	  RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	elseif fb ==  1
	   Acell(ig) =  Acell(ig) - dz*ri;
	  RAcell(ig) = RAcell(ig) - dz*ri^2/2;
	  ARcell(ig) = ARcell(ig) - dz*log(ri);
	end
      end
    elseif fa == -1 % contour went out through bottom of the cell
      if fb ~= 1 % contour didn't come in through bottom of the cell
	% Lower outer corner must be inside polygon
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if  ig-1+nz > 0 && ig+nz < ngg && ncell(ig-1+nz) == 0
	  Acell(ig-1+nz) = Ag;
	end
	if fb == -nz % contour came in through outer edge of the cell
           Acell(ig) =  Acell(ig) + (zb-zd)*ro;
          RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	else
           Acell(ig) =  Acell(ig) + dz*ro;
          RAcell(ig) = RAcell(ig) + dz*ro^2/2;
          ARcell(ig) = ARcell(ig) + dz*log(ro);
	  if fb == nz % contour came in from inner edge
             Acell(ig) =  Acell(ig) + (zb-zu)*ri;
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	  end
	end
      end
    elseif fa == 1 % contour went out through top of the cell
      if fb ~= -1 % contour didn't come in through top of the cell
	% Upper inner corner must be inside polygon
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ig+1-nz > 0 && ig+2-nz < ngg && ncell(ig+1-nz) == 0
	  Acell(ig+1-nz) = Ag;
	end
	if fb == nz % contour came in from inner edge
           Acell(ig) =  Acell(ig) + (zb-zu)*ri;
          RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	else
           Acell(ig) =  Acell(ig) - dz*ri;
          RAcell(ig) = RAcell(ig) - dz*ri^2/2;
          ARcell(ig) = ARcell(ig) - dz*log(ri);
	  if fb == -nz % contour came in through outer edge of the cell
             Acell(ig) =  Acell(ig) + (zb-zd)*ro;
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	  end
	end
      end
    end    
  end
end

% Find the rest of the cells that must be inside polygon
% and also find approximate solutions for ARcell of partly covered cells
for j = 1:nz-1
  for k = 1:nr-1
    if ncell(j,k) == 0 & Acell(j,k) > 0 % if true then j,k is inside polygon
      if ncell(j,k+1) == 0 % if true then j,k+1 is also inside polygon
	Acell(j,k+1) = Ag;
      end
      if ncell(j+1,k) == 0 % if true then j+1,k is also inside polygon
	Acell(j+1,k) = Ag;
      end
    end
  end
end
RAcell(Acell==Ag) = RA(Acell==Ag);
ARcell(Acell==Ag) = AR(Acell==Ag);

Acell(Acell < 0) = Acell(Acell < 0)+Ag;
ARcell(ARcell < 0) = ARcell(ARcell < 0)+AR(ARcell < 0);
RAcell(RAcell < 0) = RAcell(RAcell < 0)+RA(RAcell < 0);
