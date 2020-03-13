function [x2,y2,s] = interparc(x,y,n,mergeit)

if ~exist('mergeit','var')
    mergeit = 0;
end

% append extra lengths to x,y to allow for endpoint calcs
if size(x,2) ~= 1, x = x'; end
if size(y,2) ~= 1, y = y'; end
x = [x; 1];
y = [y; 1];

% arc length along the given points
arclens = cumsum([0; sqrt(diff(x).^2 + diff(y).^2)]);

% evenly distributed arc lengths
s = linspace(0,arclens(end-1),n);


for i = 1:n
    
    k = find(arclens > s(i), 1) - 1;
    dk = (s(i)-arclens(k)) / (arclens(k+1)-arclens(k)); % remainder
    
    
    x2(i) = x(k) + dk * (x(k+1) - x(k));
    y2(i) = y(k) + dk * (y(k+1) - y(k));
    
end


% Also include original x,y in outputs x2,y2
if mergeit
    v = [x2' y2' s'; x y arclens];
    v2 = sortrows(v,3);
    
    x2 = v2(:,1);
    y2 = v2(:,2);
    s = v2(:,3);            
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    