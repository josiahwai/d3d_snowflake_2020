root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,rxP,rxS,zxP,zxS] = unpack(dum.args);
xp = [rxP rxS zxP zxS];

% run the sims
try 
%   eq = designeq_ml(xp,shot,time_ms);
%   save('eq','eq')
  load('eq')
  sim = heatsim_ml(eq,shot,time_ms,1);  
  J = measure_cost(sim);
catch
  J = 10;
  fprintf('Warning: simulation did not run')
end
  
fprintf(['\nCost: ' num2str(J) '  XP: ' num2str(xp) '\n'])

% save('sim','sim');


























