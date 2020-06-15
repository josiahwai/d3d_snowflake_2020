% ========
% SETTINGS
% ========
clear
saveit = 1;
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';

% =================================
% LOAD JPAR AND XPS
% ================================

d = dir([topdir '/**/*eqs*']);

j0 = [];
jf = [];
dj = [];
dxp = [];
snowtype = {};
i_snowtpe = [];
shots = [];
times = [];
k = 0;


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
    end
  end
  catch 
  end
end


jf = real(jf);
j0 = real(j0);
dj = real(dj);

if saveit
  
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
  
%   djmax = max(abs(dj(:,52:end)'));
  
  savedir = '/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_jboot_and_xp/';
  save(  [savedir 'sim6_14'], 'sim');
end








