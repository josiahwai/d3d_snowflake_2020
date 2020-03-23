% Return the polynomial terms in the taylor expansion
% f(x,v) = l1*x + l2*v + q1*x^2 + q2*x*v
%        = c*p 
% where c is the coefficients and p = [x v x^2 x*v ... x*v^2 v^3];
         
% Reference:  D.D. Ryutov, "Local properties of the magnetic field
% in a snowflake divertor" 

function p = polyterms(x,v)

n = length(x);
p = zeros(n,10);
for k = 1:n
  p(k,:) = [1 x(k) v(k) x(k)^2 x(k)*v(k) v(k)^2 x(k)^3 ...
    x(k)^2*v(k) x(k)*v(k)^2 v(k)^3];
end


