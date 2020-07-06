% ========
% SETTINGS
% ========
clear
saveit = 1;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';
root = '/u/jwai/d3d_snowflake_2020/current';


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


for i = 1:length(d)
  try
    
  load([d(i).folder '/eqs.mat'])  % load eqs
  load([d(i).folder '/xps.mat'])  % load eqs
  
  % Load CAKE data =================
  shot = eqs{1}.shotnum;
  time = eqs{1}.time * 1000;
  load(['eq' num2str(shot) '_' num2str(time) '.mat']);
  cake_j = eq.jpar;
  
  
  % ================================
  
  valid_sim = 1;
  if length(eqs) < 3  || any(isnan(xps{end})) || ...
      any(isnan(eqs{end}.jpar))  || ~isreal(xps{end}) ||...
      ~isreal(eqs{end}.jpar) || ~isreal(eqs{2}.jpar)
    
    valid_sim = 0;
  end
  
  
  if valid_sim
    if ~any(isnan(xps{end})) && ~any(isnan(eqs{end}.jpar))
      
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
      
      shots(k) = eqs{1}.shotnum;
      times(k) = eqs{1}.time * 1000;
      j0(k,:) = eqs{2}.jpar;
      jf(k,:) = eqs{end}.jpar;
      dj(k,:) = (eqs{end}.jpar - eqs{2}.jpar);
      dxp(k,:) = xps{end} - xps{1};
      presf(k,:) = eqs{end}.pres;
      pres0(k,:) = eqs{2}.pres;
      dpres(k,:) = (eqs{end}.pres - eqs{2}.pres);
      cake_j0(k,:) = cake_j;
      
    end
  end
  catch 
  end
end


jf = real(jf);
j0 = real(j0);
dj = real(dj);



% write to struct
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
sim.cake_j0 = cake_j0;


if saveit
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_jboot_and_xp/';
  save(  [savedir 'sim_data3'], 'sim');
end








