shotlist = 155328:2:155354;

% submit batch jobs every x hours
for shot = shotlist
  try
    
    sim_sfm = 0;
    sim_sfp = 1;
    constrain_sp = 0;
    
    submit_sfmodel_jobs(shot, sim_sfm, sim_sfp, constrain_sp)
    pause(3600); % 1 hour


  catch
  end
end







