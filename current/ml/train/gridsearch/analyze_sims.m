% ---------
% SETTINGS
% ---------
close all; 

shot = 155355;
time_ms = 3900;
plotit = 0;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% --------------
% Evaluate sims
% -------------

job_topdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/'];
d = dir(job_topdir);


nJobs = length(d) - 2;
sims_that_ran = ones(nJobs,1);
J = inf(nJobs,1);
rx = zeros(nJobs,2);
zx = zeros(nJobs,2);

% b = [452   623   481   477   628   602   577   603  ];
% for k = b
for k = 1:nJobs
  jobdir = [job_topdir num2str(k) '/'];  
        
  sim_fn = [jobdir 'sim.mat'];
  d = dir([jobdir '*.out']);
  out_fn = [jobdir d.name];  
  
  if isfile(sim_fn)      % was a sim file saved?     
    load(sim_fn)
    J(k) = measure_cost4(sim);      
    rx(k,:) = sim.rx;
    zx(k,:) = sim.zx;
  elseif isfile(out_fn)  % maybe the sim ran, but no sim file was saved        
    % rmdir(jobdir,'s');
    sims_that_ran(k) = 0;
    
    fid = fopen(out_fn);
    tline = fgetl(fid);    
    cost = [];
    
    while ischar(tline)   % read outfile to see if sim ran
      if contains(tline,'Cost:')
        sims_that_ran(k) = 1;
        c = textscan(tline, '%s %f %s %f %f %f %f');
        J(k) = inf; % c{2};        
      end      
      tline = my_fgetl(fid);
    end
    fclose(fid); 
  end
  
  
  % plot heat flux  
  if plotit && isfile(sim_fn)
    struct_to_ws(sim);
        
    figure(2); clf; hold on
    set(gcf,'position',[835 591 446 282])        
    
    q_scale = 1/max([qI; qX; qO]);
    qir_scale = 1/max(qir);
    clear ylim; ylim([0 1.1])        
    
    plot(sI, qI*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)
    plot(sO, qO*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)
    plot(sX, qX*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)    
    plot(sir, qir*qir_scale, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
    
    title([num2str(k) ': ' num2str(shot) ': ' num2str(time_ms) ' ms: J = ' num2str(J(k))])
    xlabel('s [cm]')
    ylabel('Heat Flux [normalized]')
    
%     figure(3); clf; hold on
%     rx = [1.2007    1.1276];
%     zx = [-1.1932   -1.3739];
%     scatter(rx,zx,100,'k','filled')
%     scatter(sim.rx,sim.zx,'r','filled')
%     axis([1 1.3 -1.4 -1])
%     set(gcf,'position',[846 176 387 315])
  end 
end
  


% cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
% cake_eq = read_eq(shot, time_ms/1000, cake_dir);
% plot_eq(cake_eq)
% axis equal
% axis([1.0 1.5 -1.4 -0.9])
% 
% 
% [Jk,k] = mink(J,10)
% Jk = Jk/max(Jk);
% for i = length(k):-1:1
% %   scatter(rx(k(i),:),zx(k(i),:),[],Jk(i)*[1 1 1],'filled')
%   scatter(rx(k(i),:),zx(k(i),:),'r','filled')
% end


  
  
  
  
  
  
  
  






















