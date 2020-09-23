% Estimate lambda_q from how power fraction changes
%
% This is to address one of the points raised in review that the effective
% lambda_q for power splitting might differ significantly from lambda_q
% obtained from Eich fit to SP4

ccc

savedir = '/u/jwai/d3d_snowflake_2020/current/review_comments/lambda_q_estimation/';

load('all_shots_times.mat');

for i= 1:length(all_shots_times)
  fprintf( [ num2str(i) ' of ' num2str(length(all_shots_times)) '\n'])
  
  shot = all_shots_times(i,1);
  time_ms = all_shots_times(i,2);
  try
    [P2(i), P4(i), drsplit(i)] = measure_sp_power( shot, time_ms);
  catch
  end
end

save('P2', 'P2')
save('P4', 'P4')
save('drsplit', 'drsplit')

scatter(drsplit, P2 ./ (P2 + P4))





function [P2, P4, drsplit] = measure_sp_power(shot, time_ms, plotit)

if ~exist('plotit', 'var'), plotit = 0; end


efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

% simulate the efit01
eq = read_eq(shot, time_ms/1000, efit_dir);
snow = analyzeSnowflake(eq);

struct_to_ws(eq.gdata);
[rzbbbs] = contourc(rg, zg, psizr, [psibry psibry]);
[romp,k] = max(rzbbbs(1,:));
zomp = rzbbbs(2,k);

dpsi = inf;
while abs(dpsi) > .0002
  [psi, psi_r] = bicubicHermite(rg, zg, psizr, romp, zomp);
  dpsi = psi - snow.psixPL;
  romp = romp - .01 * psi_r * dpsi;
end

rsplit = romp;
dpsi = inf;
while abs(dpsi) > .0002
  [psi, psi_r] = bicubicHermite(rg, zg, psizr, rsplit, zomp);
  dpsi = psi - snow.psixSL;
  rsplit = rsplit - .01 * psi_r * dpsi;
end

drsplit = rsplit - romp;

if plotit
  plot_eq(eq)
  scatter(eq.gdata.rbbbs, eq.gdata.zbbbs)
end

% Load heat flux data q(s,t), s=distance along limiter, and t=time
root = '/u/jwai/d3d_snowflake_2020/current/';
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
[~,k] = min(abs(t-time_ms));
qperp = qperp(k,:);
slim = s;

if plotit
  figure
  plot(slim,qperp)
end


load('d3d_obj_mks_struct_6565.mat');
limdata = tok_data_struct.limdata;
slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
[rlim, zlim] = calcLimDistanceInv(slimtot - slim/100, limdata);


r2 = [1.02 1.15];
r4 = [1.15 1.34];

if plotit
  figure
  plot(rlim,qperp)
  xline(r2(1), 'r', 'linewidth', 2);
  xline(r2(2), 'r', 'linewidth', 2);
  xline(r4(2), 'r', 'linewidth', 2);
end

% subtract background heat flux
q = qperp - median(qperp);

i2 = find( r2(1) < rlim & rlim < r2(2));
i4 = find( r4(1) < rlim & rlim < r4(2));

P2 = 2 * pi * trapz( rlim(i2), rlim(i2) .* q(i2));
P4 = 2 * pi * trapz( rlim(i4), rlim(i4) .* q(i4));
end


































