% ---------
% SETTINGS
% ---------
shot = 155355;
time_ms = 3000;
plotit = 1;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));


% --------------
% Evaluate sims
% -------------

job_topdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/'];
d = dir(job_topdir);

J = [];
for k = 1:length(d)-2
  jobdir = [job_topdir num2str(k)];  
  load([jobdir '/sim.mat'])
  
  J(k) = measure_cost(sim);
  
  % plot heat flux  
  if plotit    
    struct_to_ws(sim);
        
    figure(2); clf; hold on
    
    xlabel('s [cm]')
    ylabel('Heat Flux [normalized]')
    title([num2str(shot) ': ' num2str(time_ms) ' ms: J = ' num2str(J(k))])
    
    q_scale = 1/max([qI; qX; qO]);
    qir_scale = 1/max(qir);
    ylim([0 1.1])        
    
    plot(sI, qI*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)
    plot(sO, qO*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)
    plot(sX, qX*q_scale, '-og', 'LineWidth', 1, 'MarkerSize', 2)    
    plot(sir, qir*qir_scale, '-ok', 'LineWidth', 1, 'MarkerSize', 2)
  end

  
end
  

























