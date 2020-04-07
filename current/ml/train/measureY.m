% -----------
% SETTINGS
% -----------
clear; clc; close all;

shot = 155355;
plotit = 1;
saveit = 1;

% -----------
% SCAN OUTFILE
% -----------
root = '/u/jwai/d3d_snowflake_2020/current/';
shotdir = [root 'ml/train/jobs/' num2str(shot) '/'];

times = [];
d = dir(shotdir);
for k = 1:length(d)-2
  times(k) = str2double(d(k+2).name);
end

x_all=[]; y_all=[];
for iTime = 1:length(times)
  
  time_ms = times(iTime);
  
  
  % scan outfile to find best cost
  jobdir = [shotdir num2str(time_ms) '/1/'];
  
  d = dir([jobdir '*.out']);
  fn = d.name;
  fid = fopen([jobdir fn]);
  tline = fgetl(fid);
  
  cost = []; xp = []; k = 0;
  while ischar(tline)
    if contains(tline,'Cost:')
      c = textscan(tline, '%s %f %s %f %f %f %f');
      k = k+1;
      [~,cost(k),~,rxP,rxS,zxP,zxS] = unpack(c);
      xp(k,:) = [rxP rxS zxP zxS];
      
    end
    tline = fgetl(fid);
  end
  
  cost_threshold = .002;
  [J,kBest] = min(cost);
  if J > cost_threshold
    warning('Cost appears too high')
    kBest = [];
  end
  
  xp_best = xp(kBest,:);
  
  % ---------------------
  % SIMULATE WITH BEST XP
  % ---------------------
  
  % analyze previous eq (efit)
  efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
  eq0 = read_eq(shot, time_ms/1000, efit_dir);
  sim0 = heatsim_ml(eq0.gdata,shot,time_ms,plotit);
  snow0 = analyzeSnowflake(eq0);
  
  % obtain / analyze new eq
  eq1 = designeq_ml(xp_best,shot,time_ms);
  sim1 = heatsim_ml(eq1,shot,time_ms,plotit);
  snow1 = analyzeSnowflake(eq1);
  
  if plotit
    plot_eq(eq0)
    hold on
    axis equal
    axis([1.0 1.5 -1.4 -0.9])
    plot(xp(1,1:2), xp(1,3:4),'xb', 'linewidth', 3);
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


% save data
if saveit
  savedir = [root 'ml/train/train_data/'];
  fnx = [num2str(shot) '_x.mat'];
  fny = [num2str(shot) '_y.mat'];
  
  save([savedir fnx], 'x_all');
  save([savedir fny], 'y_all');
end
  






















