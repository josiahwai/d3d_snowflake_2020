clear
shot = 155328;

load(['sol_params_et_' num2str(shot) '.mat'])
lambdaq = sol_params.coeff(:,2);


load(['sol_params_' num2str(shot) '.mat'])
lambdaq = [lambdaq; sol_params.coeff(:,2)];




plot(lambdaq)











