
ccc

iCoil = 20;
shot = 155478;


root = '/u/jwai/d3d_snowflake_2020/current/';
coil_dir = [root 'inputs/coil_data/' num2str(shot)];
coil_data = ['/coil_currents_' num2str(shot) '.mat']; % loads currents, t

load([coil_dir coil_data]);

% smooth coil currents (very noisy w.r.t time)
cc_smooth = zeros(size(coil_currents));
nc = min(size(coil_currents));

for i = 1:nc
    cc_smooth(i,:) = smooth(coil_currents(i,:), 10);
end


figure(12)
hold on
plot(t, coil_currents(iCoil,:), 'b', 'linewidth', 1);
plot(t, cc_smooth(iCoil,:), 'r', 'linewidth', 1s);
xlim([0 5000])
title(names(iCoil,:))
