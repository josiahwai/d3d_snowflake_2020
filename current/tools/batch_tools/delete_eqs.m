% ========
% SETTINGS
% ========
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/';
topdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfp_constrain_sp/';

% =====================
% ONLY SAVE THE LAST EQ
% =====================

d = dir([topdir '/**/*eqs*']);

for k = 1:length(d)
  try
    load([d(k).folder '/eqs.mat'])
    
    % eqs = eqs([1 2 end]);
    for i = 2:3
      try
      eqs{i} = rmfield(eqs{i}, {'r', 'b', 'p'});
      catch
      end
    end
    
    save([d(k).folder '/eqs.mat'], 'eqs');
  catch
    warning(d(k).folder(end-15:end))
  end
end