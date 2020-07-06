function rep = pdipmqpneq2(H,h,E,f,kmax,epstop,eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   rep = pdipmqpneq2(H,h,E,f,kmax,epstop,eta)
%
%  PURPOSE: Solve the particular quadratic program (no equality constraint)
%           min_z (0.5*x'*H*x+h'*x), s.t. E*x<=f  (QP/neq)
%
%  INPUTS: H, Hessian term nx-by-nx (decision variable x is nx-by-1)
%          h, linear term nx-by-1
%          E, inequality constraint matrix nz-by-nx (nz is # constraints)
%          f, inequality constraint rhs data nz-by-1
%          kmax, maximum number of iterations, range 1 to 200
%          epstop, relative tolerance, range 1e-12 to 1e-3 (def. 1e-9)
%          eta, dampening factor, 0<eta<1 (default 0.95)
%
%          Any/all of kmax=[],epstop=[],eta=[] is allowed and sets 
%          the parameter to the default; but this generates a message
%          that warns about defaulting behavior.
%
%  OUTPUTS: rep, structure with solution and other parameters
%	        rep.descriptions details the contents of rep
%
%  METHOD: Uses a primal-dual interior point method based on the 
%          standard predictor-corrector technique.
%          The initial point is automatically set.
%          All constraints must be present.
%          This code uses a relative norm stop condition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Erik Olofsson ON 2015-08-20
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rep=struct;

if nargin<7 || (nargin==7 && isempty(eta))
    % dampening factor eta, default value
    eta = 0.95;
    disp(['*** warning (',mfilename(),'): defaulting eta=',num2str(eta)]);
end

assert(numel(eta)==1 && eta>0 && eta<1, 'eta must be a scalar > 0 and < 1');

if nargin<6 || (nargin>=6 && isempty(epstop))
    % default value for relative tolerance
    epstop = 1e-8;
    disp(['*** warning (',mfilename(),'): defaulting epstop=',num2str(epstop)]);
end

if nargin<5 || (nargin>=5 && isempty(kmax))
    % default max iterations
    kmax = 100;
    disp(['*** warning (',mfilename(),'): defaulting kmax=',num2str(kmax)]);
end

assert(numel(kmax)==1, 'kmax must be a scalar');
assert(kmax>=1 && kmax<=200, 'Allowed range for kmax is 1 to 200');

% Do some size-consistency checking
nx = size(H,1); % # primal variables
assert(size(H,2)==nx, 'H must be a square matrix');
assert(size(h,1)==nx, 'size(h,1) must equal size(H,1)');
assert(size(h,2)==1, 'size(h,2) must be 1');
nz = size(E,1); % # inequality constraints
assert(size(E,2)==nx, 'size(E,2) must equal size(H,1)');
assert(size(f,1)==nz, 'size(f,1) must equal size(E,1)');
assert(size(f,2)==1, 'size(f,2) must be 1');

% Initial point
x=zeros(nx,1);
z=ones(nz,1);
s=ones(nz,1);

if numel(epstop)==1
    epstop=ones(3,1)*epstop;
end

assert(numel(epstop)==3, 'epstop must contain 1 or 3 elements');
assert(min(epstop)>0 && max(epstop)<=1e-3, ...
  'Allowed values in epstop are from ~1e-12 to 1e-3');
thrL=(norm(h)+1)*epstop(1);
thrs=(norm(f)+1)*epstop(2);
thrmu=epstop(3);

e=ones(nz,1);
k=0;
rL=H*x+h+E'*z;
rs=s+E*x-f;
rsz=s.*z;
mu=sum(z.*s)/nz;
while (k<=kmax && (norm(rL)>=thrL || norm(rs)>=thrs || abs(mu)>=thrmu))
    r_bar=-E'*((rsz-z.*rs)./s);
    h_bar=-(rL+r_bar);
    H_bar=H+E'*diag(z./s)*E;
    L=chol(H_bar,'lower');
    dx_a=L'\(L\h_bar);
    ds_a=-rs-E*dx_a;
    dz_a=-(rsz+z.*ds_a)./s;
    alpha_a=1; idx_z=find(dz_a<0);
    if ~isempty(idx_z)
        alpha_a=min(alpha_a,min(-z(idx_z)./dz_a(idx_z)));
    end
    idx_s=find(ds_a<0);
    if ~isempty(idx_s)
        alpha_a=min(alpha_a,min(-s(idx_s)./ds_a(idx_s)));
    end
    mu_a=((z+alpha_a*dz_a)'*(s+alpha_a*ds_a))/nz;
    sigma=(mu_a/mu)^3;
    rsz=rsz+ds_a.*dz_a-sigma*mu*e;
    r_bar=-E'*((rsz-z.*rs)./s);
    h_bar=-(rL+r_bar);
    dx=L'\(L\h_bar);
    ds=-rs-E*dx;
    dz=-(rsz+z.*ds)./s;
    alpha=1; idx_z=find(dz<0);
    if ~isempty(idx_z)
        alpha=min(alpha,min(-z(idx_z)./dz(idx_z)));
    end
    idx_s=find(ds<0);
    if ~isempty(idx_s)
        alpha=min(alpha,min(-s(idx_s)./ds(idx_s)));
    end
    x=x+eta*alpha*dx;
    z=z+eta*alpha*dz;
    s=s+eta*alpha*ds;
    k=k+1;
    rL=H*x+h+E'*z;
    rs=s+E*x-f;
    rsz=s.*z;
    mu=sum(z.*s)/nz;
end

% Output
des.x = 'Solution vector';
rep.x = x;

des.fx = 'Objective value for solution vector';
rep.fx = 0.5*x'*H*x + h'*x;

des.iters = 'Number of iterations to reach solution';
rep.iters = k;

des.epstop = 'Relative tolerance stop condition';
rep.epstop = epstop;

des.maxiters = 'Maximum number of iterations that was allowed';
rep.maxiters = kmax;

des.z = 'Lagrange multiplier';
rep.z = z;

des.s = 'Slack vector';
rep.s = s;

des.abs_mu = 'Complementarity measure';
rep.abs_mu = abs(mu);

rep.descriptions = des;

end
