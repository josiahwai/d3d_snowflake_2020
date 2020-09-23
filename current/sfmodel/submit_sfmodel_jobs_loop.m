shotlist = 155328:2:155352;

% submit batch jobs every x hours
for shot = shotlist
  try
    
    sim_sfm = 1;
    sim_sfp = 0;
    constrain_sp = 0;
    
    submit_sfmodel_jobs(shot, sim_sfm, sim_sfp, constrain_sp)
    pause(3600); % 1 hour


  catch
  end
end







