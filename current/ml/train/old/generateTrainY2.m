clear all; clc; close all;
% -------
% TESTING
% -------
shot = 155355;
time_ms = 4000;
runbatch = 1;
plotit = 1;
smax = 0.12;
ds = 0.03;
dpsi = .005;


root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% ---------------------------------------------
% DOES THE IR PREDICT SNOWMINUS (3 HF PEAKS)?
% ---------------------------------------------
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


% ----------------------
% ANALYZE SNOWFLAKE
% ----------------------

% load eq
efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
efit_eq = read_eq(shot, time_ms/1000, efit_dir);
efit_snow = analyzeSnowflake(efit_eq);
struct_to_ws(efit_eq.gdata);
[psixP0, psixS0] = unpack([efit_snow.psixPL efit_snow.psixSL]);
[rxP0,rxS0,zxP0,zxS0] = unpack([efit_snow.rx efit_snow.zx]);


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

psixS = bicubicHermite(rg,zg,psizr,rxS,zxS);
psixP = bicubicHermite(rg,zg,psizr,rxP,zxP);

k1 = abs(psixP-psixP0) < dpsi;
k2 = abs(psixS-psixS0) < dpsi;
rxP=rxP(k1); zxP=zxP(k1);
rxS=rxS(k2); zxS=zxS(k2);

iCombos = combvec(1:length(rxP), 1:length(rxS))'; 
ip = iCombos(:,1); is = iCombos(:,2);

xp_list = [rxP(ip) rxS(is) zxP(ip) zxS(is)];

if plotit
  plot_eq(efit_eq)
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
  
  batchscript = [root 'ml/train/sfd_costfun.sbatch'];
  system(['dos2unix ' batchscript]);
  
  for k = 1:length(xp_list)
    
    args = [shot time_ms xp_list(k,:)];
    jobdir = [job_topdir num2str(k)];
    
    mkdir(jobdir);
    save([jobdir '/args.mat'], 'args');
    
    % copy script to jobdir
    jobscript = [root 'ml/train/sfd_costfun.m'];
    copyfile(jobscript, jobdir);
    
    % cd and submit batch job
    cd(jobdir)
    system(['sbatch ' batchscript]);
    cd(job_topdir)
  end
  
  cd([root 'ml/train/'])
end










