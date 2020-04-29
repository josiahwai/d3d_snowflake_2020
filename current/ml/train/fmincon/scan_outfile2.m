% -----------
% SETTINGS
% -----------
clear; clc; close all;

shot = 165288;
saveit = 1;

% -----------
% SCAN OUTFILE
% -----------
root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));
shotdir = [root 'ml/train/jobs/' num2str(shot) '/'];
savedir = [root 'ml/train/train_data/'];
if ~exist(savedir,'dir'), mkdir(savedir); end

times = [];
d = dir(shotdir);
for k = 1:length(d)-2
  times(k) = str2double(d(k+2).name);
end

J0=[]; Jmin=[]; xp0=[]; xp_min=[];
for iTime = 1:length(times)
  
  time_ms = times(iTime)
  
  % scan outfile to find best cost
  jobdir = [shotdir num2str(time_ms) '/1/'];
  
  d = dir([jobdir '*.out']);
  fn = d.name;
  fid = fopen([jobdir fn]);
  tline = fgetl(fid);
  
  cost = []; xp = []; k = 0; iter = 0;
  while ischar(tline) && iter < 2000
%     iter = iter + 1;
    
    if contains(tline,'Cost:')
      c = textscan(tline, '%s %f %s %f %f %f %f');
      k = k+1;
      [~,cost(k),~,rxP,rxS,zxP,zxS] = unpack(c);
      xp(k,:) = [rxP rxS zxP zxS];
    end
    
    tline = my_fgetl(fid);
  end
  fclose(fid);
  
  
  % record the best x-pt
  J0(iTime)  = cost(1);
  xp0(iTime,:) = xp(1,:);
  
  [Jmin(iTime),k] = min(cost);
  xp_min(iTime,:) = xp(k,:);
 
end

xpts = struct('shot', shot, 'times', times', 'J0', J0', 'Jmin', Jmin', ...
  'xp0', xp0, 'xp_min', xp_min);

 
% save data
if saveit
  fn = [num2str(shot) '_xpts.mat'];
  save([savedir fn], 'xpts');
end

  
  
  
  
  
  
  
  
  