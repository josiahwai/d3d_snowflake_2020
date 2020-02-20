
fid = fopen('/u/pvail/d3d_snowflake_2019/gsedge/g165286.03500', 'r');

sline = fgetl(fid);

dum = sscanf(sline(41:end), '%f', 3);

nr = dum(2);
nz = dum(3);

% Read (neqdsk,2020) xdim, zdim, rzero, rgrid(1), zmid

dum = fscanf(fid, '%f', 5);

rdim   = dum(1);
zdim   = dum(2);

rzero  = dum(3);

rgrid1 = dum(4);
zmid   = dum(5);

% Read (neqdsk,2020) rmaxis, zmaxis, ssimag, ssibry, bcentr

dum = fscanf(fid, '%f', 5);

rmaxis = dum(1);
zmaxis = dum(2);

ssimag = dum(3);
ssibry = dum(4);

bcentr = dum(5);

% Read (neqdsk,2020) cpasma, ssimag, xdum, rmaxis, xdum

dum = fscanf(fid, '%f', 5);

cpasma = dum(1);

% Read (neqdsk,2020) zmaxis, xdum, ssibry, xdum, xdum

dum = fscanf(fid, '%f', 5);

% Read (neqdsk,2020) fpol(i) where i=1:nw

fpol = fscanf(fid, '%f', nr);

% Read (neqdsk,2020) pres(i) where i=1:nw

pres = fscanf(fid, '%f', nr);

% Read (neqdsk,2020) ffprim(i) where i=1:nw

ffprim = -sign(cpasma)*fscanf(fid, '%f', nr);

% Read (neqdsk,2020) pprime(1) where i=1:nw

pprime = -sign(cpasma)*fscanf(fid, '%f', nr);

% Read (neqdsk,2020) psizr(i,j) where i=1:nw and j=1:nh

psirz = fscanf(fid, '%f', [nr nz]);

% Read (neqdsk,2020) qpsi(i) where i=1:nw

qpsi = fscanf(fid, '%f', nr);

% Read (neqdsk,2020) nbbbs, limitr

dum = fscanf(fid, '%f', 2);

nbbbs  = dum(1);
limitr = dum(2);

% Read (neqdsk,2020) rbbbs(i), zbbbs(i) where i=1:nbbbs

dum = fscanf(fid, '%f', [2 nbbbs]);

rbbbs = dum(1,:)';
zbbbs = dum(2,:)';

% Read (neqdsk,2020) xlim(i), ylim(i) where i = 1:limitr

dum = fscanf(fid, '%f', [2 limitr]);

rlim = dum(1,:)';
zlim = dum(2,:)';

% Construct the EFIT (r,z) grid
 
rg = linspace(rgrid1, rgrid1 + rdim, nr)';
zg = linspace(zmid - zdim/2, zmid + zdim/2, nz)';

[rgg,zgg] = meshgrid(rg,zg);

% Calculate flux quantities in real units

psizr = -2*pi*psirz';

psimag = -2*pi*ssimag;
psibry = -2*pi*ssibry;

% Determine which of the grid points are inside plasma boundary

idxIn = inpolygon(rgg, zgg, rbbbs, zbbbs);

% Compute the toroidal current density on the plasma grid

mu0 = 4*pi*(1e-7);

jphi = zeros(nz,nr);

psilin = linspace(psimag, psibry, length(pprime))';

pprime_rz = interp1(psilin, pprime, psizr(idxIn), 'spline');
ffprim_rz = interp1(psilin, ffprim, psizr(idxIn), 'spline');

jphi(idxIn) = rgg(idxIn).*pprime_rz + ffprim_rz./(mu0*rgg(idxIn));
