% fit lambda_q from power split
ccc
load('P2.mat')
load('P4.mat')
load('drsplit.mat')
load('all_shots_times.mat')

iuse = find(all_shots_times(:,2) > 4000 & drsplit' > -.002);
P4 = P4(iuse); 
P2 = P2(iuse);
drsplit = drsplit(iuse);

f = double(P4 ./ (abs(P2) + P4));
scatter(drsplit, f);




power_eqn = @(lambda_q, drsplit) exp(-drsplit./lambda_q);
lambda_q0 = 0.002;
lambda_q_fit = lsqcurvefit( power_eqn, lambda_q0, drsplit, f);

r = 0:.0001:.008;
f_fit = exp(-r ./ lambda_q_fit);
hold on
plot(r,f_fit, 'linewidth', 2)

ylabel('P4 / (P4 + P2)')
xlabel('Midplane Separatrix Separation [mm]')



















