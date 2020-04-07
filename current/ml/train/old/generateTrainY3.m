% plot the uncertainty zones around the x-points
% obtain sample pts to try for x-pt locations, xp_list = [rxP rxS zxP zxS]


% -------
% TESTING
% -------
clear; clc; close all;
shot = 155328;
time_ms = 4000;
runbatch = 0;
plotit = 1;
saveit = 0;
dpsi = .005;
smax = 0.12;
ds = 0.03;

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
struct_to_ws(efit_snow);
struct_to_ws(efit_eq.gdata);

[rxP,rxS] = unpack(rx);
[zxP zxS] = unpack(zx);

% ---------------------------
% OBTAIN THE UNCERTAINTY ZONES
% ---------------------------
ng = 257;
[psizr, rg, zg] = regrid(rg, zg, psizr, ng, ng);
rgg = repmat(rg',ng,1);
zgg = repmat(zg,1,ng);

% 8 pts around grid pt, CW starting from upper left corner
neighbors = [-ng-1 -1 ng-1 ng ng+1 +1 -ng+1 -ng];


% Zone around primary x-pt
% ........................

% points within dpsi of primary x-pt and nearby
iP = find(abs(psizr-psixPL) < dpsi);

dist_from_xp = sqrt(sum(([rxP zxP] - [rgg(iP) zgg(iP)]).^2,2));
iP = iP(dist_from_xp < smax);

% eliminate pts w/o at least 4 neighbors
% (i.e., identify only the large zone(s) around the x-pts, not scattered pts)
nMin = 4;
for dum = 1:length(iP)
    nNeighbors = sum(ismember(iP+neighbors,iP),2);
    iP = iP(nNeighbors >= nMin);
    if min(nNeighbors) >= nMin, break; end
end

% Zone around primary x-pt
% ........................

% points within dpsi of secondary x-pt and nearby
iS = find(abs(psizr-psixSL) < dpsi);

dist_from_xp = sqrt(sum(([rxS zxS] - [rgg(iS) zgg(iS)]).^2,2));
iS = iS(dist_from_xp < smax);

for dum = 1:length(iS)
    nNeighbors = sum(ismember(iS+neighbors,iS),2);
    iS = iS(nNeighbors >= nMin);
    if min(nNeighbors) >= nMin, break; end
end

% Find the boundary around these pts
b = boundary(rgg(iP), zgg(iP));
bP = iP(b);
b = boundary(rgg(iS), zgg(iS));
bS = iS(b);


% Points spaced equally by ds
% ...........................

% sample x-pts
ppa = 1/ds^2;
[rP,zP] = polygrid(rgg(bP), zgg(bP), ppa);
[rS,zS] = polygrid(rgg(bS), zgg(bS), ppa);


if plotit
  plot_eq(efit_eq)
  axis([1.0 1.5 -1.4 -0.9])
  scatter(rgg(iP),zgg(iP),'b')
  scatter(rgg(iS),zgg(iS),'b')
  
  plot(rgg(bP),zgg(bP), 'g', 'linewidth', 2)
  plot(rgg(bS),zgg(bS), 'g', 'linewidth', 2)
%   scatter(rP,zP,'g','filled')
end


% -------------------------------
% COMBINATIONS OF DIFFERENT X-PTS
% ------------------------------

iCombos = combvec(1:length(rP), 1:length(rS))'; 
ip = iCombos(:,1); is = iCombos(:,2);

xp_list = [rP(ip) rS(is) zP(ip) zS(is)];









































%
%
%
%
% % ========================================
% % ANALYZE THE SNOWFLAKE EQUILIBRIUM (EFIT)
% % ========================================
%
% % Load tokamak definition
% tokdir = [root 'inputs/tok_data/'];
% tokdata = 'd3d_obj_mks_struct_129129.mat';
% load([tokdir tokdata]);
%
% % Compute total length of the limiter [m]
% limdata = tok_data_struct.limdata;
% sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
%
% % interpolate for higher point density along limiter
% [rlim,zlim] = interparc(limdata(2,:), limdata(1,:), 500, true);
%
% % distance to limiter landmarks
% s45Deg1 = sLimTot - calcLimDistance(limdata(2,79), limdata(1,79), limdata);
% s45Deg2 = sLimTot - calcLimDistance(limdata(2,78), limdata(1,78), limdata);
%
% % Configure the plots
% if plotit
%     figure(11)
%     plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
%     hold on
%     axis equal
%     axis([1.0 1.5 -1.4 -0.9])
%     xlabel('R [m]')
%     ylabel('Z [m]')
%     title([int2str(shot) ': ' int2str(time_ms) ' ms'])
% end
%
% rg = eq.rg;
% zg = eq.zg;
%
% % different naming conventions for this info
% try
%     bzero = eq.bzero;   % eq generated from gfile
%     rzero = eq.rzero;
% catch
%     bzero = eq.btsurf;  % eq generated from gsdesign
%     rzero = eq.rsurf;
% end
%
% psizr  = eq.psizr;
% psimag = eq.psimag;
% psibry = eq.psibry;
%
% % finer mesh
% ng = 257;
% [psizr, rg, zg] = regrid(rg, zg, psizr, ng, ng);
% rgg = repmat(rg',ng,1);
% zgg = repmat(zg,1,ng);
%
% % Locate the snowflake in the lower divertor
% rExp   =  1.1500;
% zExp   = -1.2500;
% rhoExp =  0.1000;
%
% [rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);
%
% [rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
% [rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);
%
% if abs(psixSL - psibry) < abs(psixPL - psibry)
%     swap(psixPL, psixSL);
%     swap(rxPL, rxSL);
%     swap(zxPL, zxSL);
% end
%
% if plotit
%     contour(rg, zg, psizr, [psixPL psixPL], 'b', 'LineWidth', 2)
%     contour(rg, zg, psizr, [psixSL psixSL], 'b', 'LineWidth', 1)
%
%     plot(rxPL, zxPL, 'xb', 'Markersize', 12, 'LineWidth', 3)
%     plot(rxSL, zxSL, 'xb', 'Markersize', 12, 'LineWidth', 3)
% end
%
%
% % =====================================
% % FIND EVENLY DISTRIBUTED SAMPLE POINTS
% % =====================================
% a = polyarea(rgg(bP),zgg(bP)) + polyarea(rgg(bS),zgg(bS));
% ppa = np/a;
%
% [rinP,zinP] = polygrid(rgg(bP), zgg(bP), ppa);
% [rinS,zinS] = polygrid(rgg(bS), zgg(bS), ppa);
%
% if plotit
%     scatter(rinP,zinP,'g')
%     scatter(rinS,zinS,'g')
% end
%
% % =============================================================
% % DOES qperp FROM IR PREDICT SNOWPLUS OR SNOWMINUS (2vs3 peaks)
% % =============================================================
% irdir = [root 'inputs/qperp/' num2str(shot) '/'];
% fn_ir = ['qperp_' num2str(shot) '.mat'];
% load([irdir fn_ir]) % loads qperp, t
%
% [~,k] = min(abs(t-time_ms));
% q = qperp(k,:);
%
% iX = find(s >= 115 & s <= 145);
%
% % find region between x-pts peak
% [~,k] = findpeaks(q(iX),'NPeaks',1,'sortstr','descend','minpeakdistance',...
%     10,'minpeakheight',0.05*max(q), 'minpeakprominence', 0.05*max(q));
% ipkX = min(iX)-1+k;
% if ~isempty(k)
%     sX = s(ipkX);
% else
%     sX = 0;
% end
%
% % plot
% if plotit
%     figure(4)
%     plot(s,q)
%     hold on
%     scatter(s(ipkX),q(ipkX),'r','filled')
% end
%
% % snowflake type
% snowPlus = false;
% if isempty(ipkX), snowPlus = true; end
% snowMin = ~snowPlus;
%
%
% % =======================================
% % RESTRICT SAMPLE X-PTS BASED ON SNOWTYPE
% % =======================================
%
% psiP = bicubicHermite(rg,zg,psizr,rinP,zinP);
% psiS = bicubicHermite(rg,zg,psizr,rinS,zinS);
%
% iCombos = combvec(1:length(rinP), 1:length(rinS))';
%
% % psiDiff = psiP(iCombos(:,1)) - psiS(iCombos(:,2));
% rDiff = rinS(iCombos(:,2)) - rinP(iCombos(:,1));
%
% % re-sort according to psi difference between sample pts
% % [psiDiff,k] = sort(psiDiff);
% [rDiff,k] = sort(rDiff);
% iCombos = iCombos(k,:);
%
% % sample at least 30% of cases, all cases that have values of
% % rxP, rxS consistent with the snowtype
%
% slack = 0.04; % 4cm of slack for snowtype
% f = 0.3;
% if snowPlus
%     % k = find(psiDiff < 0);
%     k = find(rDiff - slack < 0);
%     if length(k) < f * length(iCombos)
%         k = 1:floor(f*length(iCombos));
%     end
% elseif snowMin
%     % k = find(psiDiff > 0);
%     k = find(rDiff + slack > 0);
%     if length(k) < f*length(iCombos)
%         k = floor((1-f)*length(iCombos)):length(iCombos);
%     end
% end
% iCombos = iCombos(k,:);
%
% % save
% ip = iCombos(:,1);
% is = iCombos(:,2);
% xp_list = [rinP(ip) rinS(is) zinP(ip) zinS(is)];
%
% if saveit
%     savedir = [root 'optimize/'];
%     fn = 'xp_list';
%     save([savedir fn], 'xp_list');
% end
%
% if plotit
%     figure(11)
%     scatter(rinP(ip),zinP(ip),'g')
%     scatter(rinS(is),zinS(is),'g')
%     if saveit
%         savefig('zones.fig')
%     end
% end
%
%
% % ======================
% % SET UP FOR BATCH JOBS
% % ======================
% if runbatch
%     job_topdir = [root 'optimize/batch/'];
%     output_dir = [root 'optimize/output/'];
%
%     % clear contents
%     if exist(job_topdir,'dir'), rmdir(job_topdir,'s'); end
%     if exist(output_dir,'dir'), rmdir(output_dir,'s'); end
%     mkdir(job_topdir)
%     mkdir(output_dir)
%
%
%     batchscript = [root 'optimize/design_eq_optimize.sbatch'];
%     system(['dos2unix ' batchscript]);
%
%     for k = 1:length(xp_list)
%
%         args = [k shot time_ms xp_list(k,:)];
%         jobdir = [job_topdir num2str(k)];
%
%         mkdir(jobdir);
%         save([jobdir '/args.mat'], 'args');
%
%         % copy script to jobdir
%         jobscript = [root 'optimize/design_eq_optimize.m'];
%         copyfile(jobscript, jobdir);
%
%         % cd and submit batch job
%         cd(jobdir)
%         system(['sbatch ' batchscript]);
%         cd(job_topdir)
%     end
%
%     cd([root 'optimize'])
% end





























