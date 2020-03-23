% Return the polynomial terms in the taylor expansion
% f(x,v) = l1*x + l2*v + q1*x^2 + q2*x*v

% dpsi/dx = c*px
% dpsi/dv = c*pv;
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

