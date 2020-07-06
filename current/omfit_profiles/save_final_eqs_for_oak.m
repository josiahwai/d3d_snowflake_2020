d = dir('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm');
% d = dir('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp_constrain_sp/155354_sfp_sp');
d(1:2) = [];
shot = 155354;
savedir = '/u/jwai/d3d_snowflake_2020/current/profiles/eqfiles_tf/';

for k = 1:length(d)
  
  time_ms = str2num(d(k).name(end-3:end));      
  load( [d(k).name '/eqs.mat']);
  
  if length(eqs) >= 3    
    eq = eqs{end};    
    fn = [savedir 'eq' num2str(shot) '_' num2str(time_ms)];    
    save(fn, 'eq');
  else
    warning(num2str(time_ms))
  end
  
end
  
 








