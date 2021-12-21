% fit lambda_q from power split

ccc
drsplit_thresh = 0.00;
% shot = 155328;
% load(['sim' num2str(shot) '.mat'])
load('sim_all.mat')
struct_to_ws(sim);

% load('/u/jwai/d3d_snowflake_2020/current/lambdaq_split/old/lambda_q_estimation/P2.mat')
% load('/u/jwai/d3d_snowflake_2020/current/lambdaq_split/old/lambda_q_estimation/P4.mat')
% load('/u/jwai/d3d_snowflake_2020/current/lambdaq_split/old/lambda_q_estimation/drsplit.mat')

f = double(P4 ./ (abs(P2) + P4));
f0 = f;
drsplit0 = drsplit;

iuse = find(drsplit' > drsplit_thresh);
[~,ioutl] = rmoutliers([f(iuse)' drsplit(iuse)']);
ioutl = false(size(iuse));

f = f(iuse(~ioutl));
drsplit = drsplit(iuse(~ioutl));


lambdaq0 = 0.001;
shift0 = 0;
param0 = [lambdaq0 shift0];

params_fit = lsqcurvefit(@power_eqn, param0, f, drsplit);

[lambdaq_fit, shift_fit] = unpack(params_fit)
mean(lambdaq_o)

% lambdaq_fit = .006

% plot fit
r = -.005:.0001:.005;
f_fit = exp(- (r + shift_fit) ./ lambdaq_fit);
f_fit(f_fit > 1) = 1;

scatter(drsplit0*1000, f0, 'r');
hold on
scatter(drsplit*1000, f, 'b');
plot(r*1000, f_fit, 'linewidth', 2)
axis([-3.8652    6.6064         0    1.1000])
ylabel('P4 / (P4 + P2)')
xlabel('Midplane Separatrix Separation [mm]')
set(gcf, 'position', [621 53 484 379])


% POWER EQUATION
function drsplit = power_eqn(params, f) 
  [lambdaq, shift] = unpack(params);
  drsplit = -lambdaq * log(f) - shift;
end


















