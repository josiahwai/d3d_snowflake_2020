ccc
plotit = 0;
% shotlist = [155330:2:155352];
shotlist = [155329, 155331, 155333];

for ishot = 1:length(shotlist)
  shot = shotlist(ishot);
  
  load(['q' num2str(shot) '.mat'])
  times = double(t(t>0));
  
  % Estimate lambda_q from how power fraction changes
  %
  % This is to address one of the points raised in review that the effective
  % lambda_q for power splitting might differ significantly from lambda_q
  % obtained from Eich fit to SP4
  
  savedir = '/u/jwai/d3d_snowflake_2020/current/lambdaq_split/';
  
  disp(shot)
  for i= 1:length(times)
    fprintf( [ num2str(i) ' of ' num2str(length(times)) '\n'])
    
    time_ms = times(i);
    
    try
      [P2(i), P4(i), drsplit(i), lambdaq_i(i), lambdaq_o(i)] = ...
        measure_sp_power( shot, time_ms, plotit);
    catch
      warning(num2str(time_ms))
    end
    
  end
  
  sim.P2 = P2;
  sim.P4 = P4;
  sim.drsplit = drsplit;
  sim.lambdaq_i = lambdaq_i;
  sim.lambdaq_o = lambdaq_o;
  save([savedir 'sim' num2str(shot)], 'sim')
  % scatter(drsplit, P4 ./ (P2 + P4))

end




function [P2, P4, drsplit, lambdaq_i, lambdaq_o] = measure_sp_power(shot, time_ms, plotit)

if ~exist('plotit', 'var'), plotit = 0; end

efit_dir = ['/p/omfit/users/jwai/projects/lambdaq_split/EFITtime/OUTPUTS/' num2str(shot) ...
  '/gEQDSK'];

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
load(['q' num2str(shot) '.mat'])
[~,k] = min(abs(t - time_ms));
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

ef = eich_fitter_dev(s', qperp', eq, tok_data_struct);
lambdaq_i = ef.fiti.lambdaQ / 100; % [m]
lambdaq_o = ef.fito.lambdaQ / 100; % [m]


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
q(q<0) = 0;

ishadow = find(slim > 160.5 & slim < 185);
q(ishadow) = 0;



i2 = find( r2(1) < rlim & rlim < r2(2));
i4 = find( r4(1) < rlim & rlim < r4(2));

P2 = 2 * pi * trapz( rlim(i2), rlim(i2) .* q(i2));
P4 = 2 * pi * trapz( rlim(i4), rlim(i4) .* q(i4));


end












































