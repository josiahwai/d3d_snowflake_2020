% minimize over x: 2-norm (A*x-b) subject to Aeq*x = beq
% returns x and lagrange multipliers, lam

function [x,lam] = leastsquare(A,b,Aeq,beq)

if nargin == 2
  Aeq=[]; beq=[];
end

nc = length(beq);

bhat = [A'*b; beq];
Ahat = [A'*A Aeq'; Aeq zeros(nc)];

xlam = pinv(Ahat)*bhat;

x = xlam(1:end-nc);
lam = xlam(end-nc+1:end);

