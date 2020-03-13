function mylegend(labels, ls, co,location)

if ~exist('location', 'var')
    location = 'best';
end

n = length(labels);
ha = zeros(n,1);

for k = 1:n
    ha(k) = plot(NaN,NaN, 'linestyle', ls{k}, 'color', co{k}, 'linewidth', 2);
end

l = legend(ha, labels, 'location', location);



