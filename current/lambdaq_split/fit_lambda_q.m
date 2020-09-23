% fit lambda_q from power split
ccc
saveit = 0;
drsplit_thresh = 0.00;
% shot = 155336
% load(['sim' num2str(shot) '.mat'])
load('sim_all.mat')
savedir = '/u/jwai/d3d_snowflake_2020/current/lambdaq_split/';

struct_to_ws(sim);

f = double(P4 ./ (abs(P2) + P4));
f0 = f;
drsplit0 = drsplit;

iuse = find(drsplit' > drsplit_thresh);
[~,ioutl] = rmoutliers([f(iuse)' drsplit(iuse)']);
% ioutl = false(size(f(iuse)));

f = f(iuse(~ioutl));
drsplit = drsplit(iuse(~ioutl));


power_eqn = @(params, drsplit) exp(-(drsplit + params(2)) ./ params(1));

lambdaq0 = 0.002;
shift0 = .001;
param0 = [lambdaq0 shift0];

params_fit = lsqcurvefit( power_eqn, param0, drsplit, f);

lambdaq_fit = params_fit(1)
shift_fit = params_fit(2)
mean(lambdaq_o)
mean(lambdaq_i)


r = -.005:.0001:.007;
f_fit = exp(- (r + shift_fit) ./ lambdaq_fit);
f_fit(f_fit > 1) = 1;


co = 	[0.85 0.32 0.1];
scatter( (drsplit0 + shift_fit)*1000, f0, 'markeredgecolor', co);
hold on
scatter( (drsplit + shift_fit) *1000, f, 'markeredgecolor', co);
plot( (r + shift_fit) *1000, f_fit, '--k', 'linewidth', 2)
axis([-3.9252    3.8   0.15    1.15])
yticks([0:0.2:1.2])
ylabel('\textrm{SP4 Power Fraction}, P4 / (P2 + P4)', 'fontsize', 12, 'interpreter', 'latex')
xlabel('$\textrm{Midplane Separatrix Separation}, r_{split}+r_{shift}\;\;  \textrm{[mm]} $', 'interpreter', 'latex')
title('Snowflake Minus Power Splitting')
legend('155330-155340')
set(gcf, 'position', [428 319 484 379])

if saveit
  saveas(gcf, [savedir 'lambdaq_eff.eps'], 'epsc')
end














