shot = '165288';
plotit = 1;
saveit = 1;

load([shot '/coil_currents.mat'])

cc = zeros(size(coil_currents));

for i = 1:18
    cc(i,:) = smooth(coil_currents(i,:), 80);
end

if plotit 
    plot(t/1000,coil_currents(4,:))
    hold on
    plot(t/1000, cc(4,:))
    hold off
    xlim([0 6])
end
    
if saveit
    save([shot '/t'], 't')
    save([shot '/names'], 'names')
    save([shot '/coil_currents_smooth'], 'cc')
end


