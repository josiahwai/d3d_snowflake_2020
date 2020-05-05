% ---------
% SETTINGS
% ---------
close all; 

shot = 155353;
plotit = 0;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% --------------
% Evaluate sims
% -------------

shotdir = [root 'ml/train/jobs/' num2str(shot) '/'];
shotdir_info = dir(shotdir);


nJobs = length(shotdir_info) - 2;
sims_that_ran = zeros(nJobs,1);
J = inf(nJobs,1);
rx = zeros(nJobs,2);
zx = zeros(nJobs,2);
times = [];
nlines = zeros(nJobs,1);

for k = 1:nJobs
  time_ms = str2num(shotdir_info(k+2).name);
  times = [times; time_ms];
  
  jobdir = [shotdir num2str(time_ms) '/'];  
        
  sim_fn = [jobdir 'sim.mat'];
  d = dir([jobdir '*.out']);
  out_fn = [jobdir d(end).name];  
  
  if isfile(sim_fn)      % was a sim file saved?     
    sims_that_ran(k) = 1;
    load(sim_fn)
    J(k) = measure_cost4(sim);      
  end
  
  
  if isfile(out_fn)  % maybe the sim ran, but no sim file was saved            
    
    fid = fopen(out_fn);
    tline = fgetl(fid);    
    cost = [];
    
    while ischar(tline)   % read outfile to see if sim ran
      if contains(tline,'Cost:')
        nlines(k) = nlines(k) + 1;        
        c = textscan(tline, '%s %f %s %f %f %f %f');
        J(k) = min(J(k), c{2});        
      end      
      tline = my_fgetl(fid);
    end
    fclose(fid); 
  end
  
end

sims_that_ran = boolean(sims_that_ran);

times(sims_that_ran);
times(~sims_that_ran);
nlines

  
  
  
  
  
  
  






















