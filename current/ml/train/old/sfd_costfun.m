dum = load('args');

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

[shot,time_ms,rxP,rxS,zxP,zxS] = unpack(dum.args);
xp = [rxP rxS zxP zxS];

try 

  % --------------------
  % RUN THE SIMULATIONS
  % --------------------
  eq = designeq_ml(xp,shot,time_ms);
  sim = heatsim_ml(eq,shot,time_ms);  
  J = measure_cost(sim);
  
catch
  J = 10;
  fprintf('Warning: simulation did not run')
end
  
fprintf(['\nCost: ' num2str(J) '  XP: ' num2str(xp) '\n'])

output = [J xp];
save('output', 'output');



























