function J = sfd_costfunction(xp,shot,time_ms,opts)

if nargin ~= 4, opts.savesim = 0; end

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




























