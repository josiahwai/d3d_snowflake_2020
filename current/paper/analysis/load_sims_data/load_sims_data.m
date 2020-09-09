% ========
% SETTINGS
% ========
clear
saveit = 1;
% topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp/';
% topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp/155330_sfp/';
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm_large_lambdaq';
root = '/u/jwai/d3d_snowflake_2020/current';
save_fn = 'sims_sfm_large_lambdaq';

cd(topdir)

% =================================
% LOAD JPAR AND XPS
% ================================
addpath(genpath(root));

d = dir([topdir '/**/*eqs*']);

cake_j0 = [];
j0 = [];
jf = [];
dj = [];
dxp = [];
snowtype = {};
i_snowtpe = [];
shots = [];
times = [];
pres0 = [];
presf = [];
dpres = [];
k = 0;
V0 = [];
Vf = [];
dzmaxix = [];
dtriu = [];
dtril = [];
ddpsi = [];


for i = 1:length(d)
  try
    
  load([d(i).folder '/eqs.mat'])  % load eqs
  load([d(i).folder '/xps.mat'])  % load eqs
  
  valid_sim = 1;
  if length(eqs) < 3  || any(isnan(xps{end})) || ...
      any(isnan(eqs{end}.jpar))  || ~isreal(xps{end}) ||...
      ~isreal(eqs{end}.jpar) || ~isreal(eqs{2}.jpar)
    
    valid_sim = 0;
  end
  
  
  if valid_sim
      k
      k = k+1;
     
      if contains(d(i).folder, 'sfp_sp')
        i_snowtype(k) = 1;
        snowtype{k} = 'sfp_sp';
      elseif contains(d(i).folder, 'sfp')
        i_snowtype(k) = 2;
        snowtype{k} = 'sfp';
      else
        i_snowtype(k) = 3;
        snowtype{k} = 'sfm';
      end
      
%       snow0 = analyzeSnowflake(eqs{2});
%       snowf = analyzeSnowflake(eqs{end});
%       ddpsi(k) = (snowf.psixPL - snowf.psixSL) - (snow0.psixPL - snow0.psixSL);
      
      
      
      
%       ddpsi(k) = (snowf.psixPL - snowf.psixSL) - (snow0.psixPL - snow0.psixSL);
      shots(k) = eqs{1}.shotnum;
      times(k) = eqs{1}.time * 1000;
      j0(k,:) = eqs{2}.jpar;
      jf(k,:) = eqs{end}.jpar;
      presf(k,:) = eqs{end}.pres;
      pres0(k,:) = eqs{2}.pres;
      dj(k,:) = (eqs{end}.jpar - eqs{2}.jpar);
      dxp(k,:) = xps{end} - xps{1};
      dpres(k,:) = (eqs{end}.pres - eqs{2}.pres);
      V0(k) = find_plasma_volume(eqs{2});
      Vf(k) = find_plasma_volume(eqs{end});
      dzmaxis(k) = eqs{3}.zmaxis - eqs{2}.zmaxis;
      dtriu(k) = eqs{3}.triu - eqs{2}.triu;
      dtril(k) = eqs{3}.tril - eqs{2}.tril;
  end
  catch 
  end
end


% write to struct
sim.i_snowtype = i_snowtype;
sim.snowtype = snowtype;
sim.shots = shots;
sim.times = times;
sim.i_snowtype = i_snowtype;
sim.snowtype = snowtype;
sim.j0 = j0;
sim.jf = jf;
sim.dj = dj;
sim.dxp = dxp;
sim.j0max = max(j0(:,52:end)');
sim.jfmax = max(jf(:,52:end)');
sim.djmax =  sim.jfmax - sim.j0max;
sim.psin = eqs{end}.psibar;
sim.pres0 = pres0;
sim.presf = presf;
sim.dpres = dpres;
sim.V0 = V0;
sim.Vf = Vf;
sim.dV = Vf - V0;
sim.dzmaxis = dzmaxis';
sim.dtriu = dtriu';
sim.dtril = dtril';
sim.ddpsi = ddpsi';

if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/';
  save(  [savedir save_fn], 'sim');
end








