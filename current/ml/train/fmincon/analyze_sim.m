% ---------
% SETTINGS
% ---------
close all; 

shot = 155354;
plotit = 0;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% ===================
% FIND SIMS TO PLOT
% ==================

shotdir = [root 'ml/train/jobs/' num2str(shot) '/'];
shotdir_info = dir(shotdir);


nJobs = length(shotdir_info) - 2;
sims_that_ran = zeros(nJobs,1);

J = inf(nJobs,1);
times = [];
nlines = zeros(nJobs,1);
sims = {};
eqs = {};
xsol = {};

for k = 1:nJobs
  time_ms = str2num(shotdir_info(k+2).name);
  times = [times; time_ms];
  
  jobdir = [shotdir num2str(time_ms) '/'];  
        
  sim_fn = [jobdir 'sim.mat'];
  eq_fn = [jobdir 'eq.mat'];
  d = dir([jobdir '*.out']);
  out_fn = [jobdir d(end).name];  
  
  % was a sim file saved? 
  % .....................
  if isfile(sim_fn)      
    sims_that_ran(k) = 1;
    load(sim_fn)
    load(eq_fn)
    sims{k} = sim;
    eqs{k} = eq;
    % J(k) = measure_cost4(sim);      
  end
 
  
  % maybe the sim ran, but no sim file was saved  
  % ............................................
  if isfile(out_fn)        
    fid = fopen(out_fn);
    tline = fgetl(fid);    
    cost = [];
    
    % read outfile
    while ischar(tline)   
      if contains(tline,'Cost:')
        nlines(k) = nlines(k) + 1;        
        c = textscan(tline, '%s %f %s %f %f %f %f');
        J(k) = min(J(k), c{2});  
        xsol{k} = [c{4} c{5} c{6} c{7}];
      end      
      tline = my_fgetl(fid);
    end
    fclose(fid); 
  end
  
end


% sims_that_ran = boolean(sims_that_ran);
% times(sims_that_ran)
% times(~sims_that_ran)


% ================
% ANALYZE AND PLOT
% ================

[dum,sims_to_plot] = mink(J,1);
% [dum,sims_to_plot] = maxk(nlines,1);
% sims_to_plot = find(~cellfun(@isempty,sims))

for k = sims_to_plot
  
  % simulate initial eq
  % ...................
  
  % load 
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  eq0 = read_eq(shot, time_ms/1000, cake_dir);
  eq0 = eq0.gdata;
  cake_snow = analyzeSnowflake(eq0);
  xp0 = [cake_snow.rx cake_snow.zx];
  
  % simulate
  sim0 = heatsim_ml(eq0,shot,time_ms);  
  J0 = measure_cost4(sim0);
  
  
  % simulate final soln
  % ...................
  
  % eq1 = eqs{k};
  % sim1 = sims{k};
  xp1 = xsol{k};
  eq1 = designeq_ml(xp1,shot,time_ms);
  sim1 = heatsim_ml(eq1,shot,time_ms);  
  J1 = measure_cost4(sim1);
  
  
  % Plot comparison
  % ===============
  
  % plot eq
  % .......
  figure(1)
  plot_eq(eq1)
  axis([0.95 1.45 -1.45 -0.95])
  plot(xp1(:,1), xp1(:,3), 'xr', 'Markersize', 12, 'linewidth',4)
  plot(xp1(:,2), xp1(:,4), 'xr', 'Markersize', 12, 'linewidth',4)

  %   contour(eq0.rg, eq0.zg, eq0.psizr, [cake_snow.psixPL cake_snow.psixSL], 'b', ...
  %     'linewidth', 1)
  plot(xp0(:,1), xp0(:,3), 'xb', 'Markersize', 12, 'linewidth',4)
  plot(xp0(:,2), xp0(:,4), 'xb', 'Markersize', 12, 'linewidth',4)
  
  set(gcf,'position', [427 357 527 336])
  
  % plot heat flux
  % ..............
  figure(2)
  hold on
  
  % plot final heat flux
  struct_to_ws(sim1);
  qir = qir /sum(qirmax); % renormalize so that sum(pks) = 1
  plot(sir,qir,'k','linewidth',1.5)
  
  s = [sI sX sO]; 
  q = [qI; qX; qO]; 
  q = q / sum(qmax);  % renormalize
  plot(s, q, 'r', 'linewidth', 1.5)  
  
  
  % plot initial heat flux
  struct_to_ws(sim0);
   
  s = [sI sX sO];
  q = [qI; qX; qO];
  q = q / sum(qmax);  % renormalize
  plot(s, q, 'b', 'linewidth', 1)
  
  axis([0.8 1.8 0 0.5])
  
  set(gcf,'position', [419 82 535 190])
end
  
  
  
  






















