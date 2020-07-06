function [xderiv,tderiv]=deriv(x,t,method)
 %
%  SYNTAX:  [xderiv,tderiv] = deriv(x,t,method)
%
%  PURPOSE:  Calculate approximate derivative of vector or matrix.
%
%  INPUT:
%    x      =  column vector or matrix to differentiate
%    t      =  times corresponding to x data (column vector), or
%               time between samples (scalar)
%    method = differentiation method (optional, default=1)
%             1 = forward difference
%             2 = backward difference
%             3 = average of forward and backward difference
%             4 = produces correct derivative for times BETWEEN input samples
%               (Strictly speaking, derivative is undefined at sample points for 
%                   piecewise linear signal)
%
%  OUTPUT:
%    xderiv = approximate derivative vector = dx/dt, at input times if method<4,
%			at constructed time vector tderiv if method=4.
%    tderiv = input t if method<4.  Strictly contains t if method=4.
%
%  If x is a matrix, each column must represent data sampled at times given
%  by vector t.
%
%  RESTRICTIONS: 
%    (1) Time t must column vector with values in ascending order.

%    (2) Size of time must be one of the following:
%	(a) same size matrix as x, i.e. columns correspond
%	(b) column vector of length = number of rows in x 
%		(time vector used for all columns of x)
%	(c) scalar = time between rows of x

%  METHOD:  
%
%  WRITTEN BY:  Mike Walker     ON      1/8/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)deriv.m	1.4 10/13/10

% TO DO:
%	- add logic to allow t to be a matrix or scalar delta-t

% check arguments OK
sx = size(x);
if(sx==0)
   disp('ERROR deriv: size of input data is zero')
   xderiv=[];
   return;
end

st = size(t);
if(st(2) > st(1))
   fprintf(' time must be column vector or scalar - no calculation\n')
   return
end

% if time is a vector, length must match data stored in x.
if any(st~=[1 1]) & st(1)~=sx(1)
   fprintf('deriv: number of rows in x must match length(t) - no calculation\n')
   return
end
if nargin < 3, method=1;, end;

if any(st~=[1 1]) 	% if t is a vector, dt must be calculated for each step
   temp = diff(t);
   dt = min(temp);
   for k=1:sx(2)
      xderiv(:,k) = diff(x(:,k))./temp;
   end
   clear temp
else
   xderiv = diff(x)./t;		% if t not a vector, t must be = dt
   dt = t;
   t = t*[0:1:size(xderiv,1)];
end

if method==1
   tderiv = t;
   xderiv(sx(1),:) = xderiv(sx(1)-1,:);
elseif method==2
   tderiv = t;
   xderiv(2:sx(1),:) = xderiv(1:sx(1)-1,:);
   xderiv(1,:) = xderiv(2,:);
elseif method==3
   tderiv = t;
   tempb = xderiv(2,:);
   tempe = xderiv(sx(1)-1,:);
   xderiv(2:sx(1)-1,:) = (xderiv(1:sx(1)-2,:)+xderiv(2:sx(1)-1,:))/2;
   xderiv(1,:) = tempb;
   xderiv(sx(1),:) = tempe;   
elseif method==4
   dt = 1e-3*dt;
   temp = xderiv;
   nsamp = size(t,1);
   tderiv = zeros(2,1);		% force it to be column vector
   for k=1:nsamp-1
      xderiv(2*k-1)=temp(k);
      xderiv(2*k)=temp(k);
      tderiv(2*k-1)=t(k);
      tderiv(2*k)=t(k+1)-dt;
   end
   tderiv(2*nsamp-2) = t(nsamp);
end

return
