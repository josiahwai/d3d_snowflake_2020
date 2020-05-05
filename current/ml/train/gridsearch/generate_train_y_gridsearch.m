clear all; clc; close all; warning('off','all')

% ---------
% SETTINGS
% ---------
shot = 155355;
time_ms = 3900;
runbatch = 0;
plotit = 1;
smax = 0.1; % grid radius
ds = 0.03;   % grid spacing

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% ------------------------------
% DEFINE GRID NEAR CURRENT X-PTS
% -------------------------------

% load eq
cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
cake_eq = read_eq(shot, time_ms/1000, cake_dir);
cake_snow = analyzeSnowflake(cake_eq);
struct_to_ws(cake_eq.gdata);
[psixpb, psixS0] = unpack([cake_snow.psixPL cake_snow.psixSL]);
[rxP0,rxS0,zxP0,zxS0] = unpack([cake_snow.rx cake_snow.zx]);

% circle around x-pts
r=[]; z=[]; k=0;
for th = linspace(0,2*pi,20)
  k = k+1;
  r(k) = smax*cos(th);
  z(k) = smax*sin(th);
end
  
% sample x-pts
ppa = 1/ds^2;
[rxP,zxP] = polygrid(rxP0+r, zxP0+z,ppa);
[rxS,zxS] = polygrid(rxS0+r, zxS0+z,ppa);


% eliminate some points that are further away in psi
psixS = bicubicHermite(rg,zg,psizr,rxS,zxS);
psixP = bicubicHermite(rg,zg,psizr,rxP,zxP);

dpsi = mean(abs([psixP - psixpb; psixS-psixS0]));

k1 = abs(psixP-psixpb) < 1.5*dpsi;
k2 = abs(psixS-psixS0) < 1.5*dpsi;
rxP=rxP(k1); 
zxP=zxP(k1);
rxS=rxS(k2); 
zxS=zxS(k2);


% combinations of the 2 possible x-pts
iCombos = combvec(1:length(rxP), 1:length(rxS))'; 
ip = iCombos(:,1); 
is = iCombos(:,2);

xp_list = [rxP(ip) rxS(is) zxP(ip) zxS(is)];


% eliminate x-pt combinations that don't match predicted snowflake type 
% .....................................................................

% Does the IR predict snowminus (3 heat flux peaks?)
irdir = [root 'inputs/qperp/' num2str(shot) '/'];
fn_ir = ['qperp_' num2str(shot) '.mat'];
load([irdir fn_ir]) % loads qperp, t

[~,k] = min(abs(t-time_ms));
q = qperp(k,:);

iX = find(s >= 115 & s <= 145);

% is there a middle heat flux peak 
[~,k] = findpeaks(q(iX),'NPeaks',1,'sortstr','descend','minpeakdistance',...
    10,'minpeakheight',0.05*max(q), 'minpeakprominence', 0.05*max(q));

snowPlus = false;
if isempty(k), snowPlus = true; end
snowMin = ~snowPlus;

iUse = [];
th = [];
for k = 1:length(xp_list)
  % check angle between (rmaxis,zmaxis), (rxP,zxP), and (rxS,zxS)
  % large angles (e.g., > 135 deg) should be snowplus
  
  pa = [cake_eq.gdata.rmaxis cake_eq.gdata.zmaxis]; 
  pb = xp_list(k,[1 3]); % primary 
  pc = xp_list(k,[2 4]); % secondary
  
  nab = (pa - pb) / norm(pa - pb);  % Normalized vectors
  nbc = (pc - pb) / norm(pc - pb);
  
  th_abc = atan2(norm(det([nab; nbc])), dot(nab, nbc)) * 180/pi;

  th(k) = th_abc;
  
%   P = [pa; pb; pc];
%   plot(P(:,1), P(:,2))
%   hold on
%   scatter(P(:,1), P(:,2),'filled')

  if snowPlus && th_abc > 135
    iUse = [iUse k];
  elseif snowMin && th_abc < 150
    iUse = [iUse k];
  end
end

xp_list = [rxP0,rxS0,zxP0,zxS0; xp_list(iUse,:)];

disp(['There are ' num2str(length(xp_list)) ' x-pt combos.'])

if plotit
  plot_eq(cake_eq)
  axis([1.0 1.5 -1.4 -0.9])
  scatter(rxS,zxS,'filled')
  scatter(rxP,zxP,'filled')
end



% ====================
% RUN BATCH JOBS
% =====================
if runbatch
  
  % clear contents
  job_topdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/'];
  if exist(job_topdir,'dir'), rmdir(job_topdir,'s'); end
  mkdir(job_topdir)
  
  batchscript = [root 'ml/train/gridsearch/sim_gridsearch.sbatch'];
  
  for k = 1:length(xp_list)
    
    args = [shot time_ms xp_list(k,:)];
    jobdir = [job_topdir num2str(k)];
    
    mkdir(jobdir);
    save([jobdir '/args.mat'], 'args');
    
    % copy script to jobdir
    jobscript = [root 'ml/train/gridsearch/sim_gridsearch.m'];
    copyfile(jobscript, jobdir);
    
    % cd and submit batch job
    cd(jobdir)
    system(['sbatch ' batchscript]);
    cd(job_topdir)
  end
  
  cd([root 'ml/train/'])
end










