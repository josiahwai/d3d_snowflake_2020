function J = sfd_costfunction(xp,shot,time_ms,saveit)

if ~exist('saveit','var'), saveit = 0; end

try 

  % --------------------
  % RUN THE SIMULATIONS
  % --------------------
  eq = designeq_ml(xp,shot,time_ms);
  sim = heatsim_ml(eq,shot,time_ms);  
  J = measure_cost4(sim);
  
  if saveit
    save('eq','eq');
    save('sim','sim');
  end
  
catch
  J = 100;
  fprintf('Warning: simulation did not run')
end
  
fprintf(['\nCost: ' num2str(J) '  XP: ' num2str(xp) '\n'])




























