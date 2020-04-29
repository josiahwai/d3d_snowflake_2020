% -----------
% SETTINGS
% -----------
clear; clc; close all;

shot = 155355;
plotit = 1;
saveit = 0;

% -----------
% SCAN OUTFILE
% -----------
root = '/u/jwai/d3d_snowflake_2020/current/';
xpts_fn = [root 'ml/train/train_data/' num2str(shot) '_xpts.mat'];

load(xpts_fn)
x_all = [];
y_all = [];
for k = 4 % 1:length(xpts.times)    
  
  if xpts.Jmin(k) < .004  % only include good sims
  
    % ---------------------
    % SIMULATE WITH BEST XP
    % ---------------------
    time_ms = xpts.times(k);
    xp0 = xpts.xp0(k,:);
    xp_best = xpts.xp_min(k,:);
    
    % analyze previous eq (efit)
    cake_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
    eq0 = read_eq(shot, time_ms/1000, cake_dir);
    sim0 = heatsim_ml(eq0.gdata,shot,time_ms,plotit);
    snow0 = analyzeSnowflake(eq0);
    
    % obtain / analyze new eq
    eq1 = designeq_ml(xp_best,shot,time_ms);
    sim1 = heatsim_ml(eq1,shot,time_ms,plotit);
    snow1 = analyzeSnowflake(eq1);
    
    if plotit
      figure(9); clf; hold on      
      plot_eq(eq0)
      plot_eq(eq1)
      axis equal
      axis([1.0 1.5 -1.4 -0.9])
      plot(xp0(1,1:2), xp0(1,3:4),'xr', 'linewidth', 3);
      plot(xp_best(1:2), xp_best(3:4),'xb', 'linewidth', 3);
    end
    
    % ---------------------
    % DEFINE TARGET VECTOR Y
    % ---------------------
    z0 = [snow0.rx snow0.zx snow0.rSnow snow0.zSnow snow0.drSnow ...
      snow0.dzSnow  snow0.sSP_heat];
    z1 = [snow1.rx snow1.zx snow1.rSnow snow1.zSnow snow1.drSnow ...
      snow1.dzSnow snow1.sSP_heat];
    
    dz = z1-z0;  % may be better at predicting changes (dz) than z itself
    y = [z1 dz];
    
    y(isnan(y))=0;
    
    % --------------------------
    % DEFINE PREDICTOR VECTOR X
    % --------------------------
    
    % features of the eq (z0) and of the heat sim
    x = [z0 sim0.s_qmax sim0.s_qirmax sim0.qmaxN sim0.qirmaxN ...
      sim0.Apk_qmaxN sim0.Apk_qirmaxN];
    
    x(isnan(x)) = 0;
    
    x_all = [x_all x'];
    y_all = [y_all y'];
  end
end


% save data
if saveit
  savedir = [root 'ml/train/train_data/'];
  fnx = [num2str(shot) '_x.mat'];
  fny = [num2str(shot) '_y.mat'];
  
  save([savedir fnx], 'x_all');
  save([savedir fny], 'y_all');
end























