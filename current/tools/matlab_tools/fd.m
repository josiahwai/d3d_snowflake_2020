% First order finite differences for vectors

function [dydx] = fd(y,x)
    
    dx = x(2) - x(1);
    n = length(x);

    if size(y,2) ~= 1, y = y'; end
    if size(x,2) ~= 1, x = x'; end
        
    
    dydx(1) = (y(2) -   y(1)) / dx;
    dydx(n) = (y(n) - y(n-1)) / dx;
    
    for i = 2:n-1
        dydx(i) = (y(i+1) - y(i-1)) / (2*dx);
    end

    
end