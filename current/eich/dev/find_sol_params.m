% Finds the thermal diffusivity chi and scrape off layer power
% decay lengths lambdaq for the inboard and outboard divertor. Takes an 
% average over the whole shot. Params are found from the eich profile

% lambdaq in [cm], chi in Wb^2/s

function [chi, lambdaq_i, lambdaq_o] = find_sol_params(shot)

coeff = [];
load('d3d_obj_mks_struct_6565.mat')

% load heat flux
root = '/u/jwai/d3d_snowflake_2020/current/';
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t



for iTime = 1:length(t)
  try
    
    time_ms = t(iTime);
    q = qperp(iTime,:)';
    
    % load eq and calculate flux expansion
    eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
    eq = read_eq(shot, time_ms, eqdir);
    
    
    [rsp, zsp, ssp, chi, fiti, fito] = eich_fitter(s, q, eq, tok_data_struct);
    
    coeff = [coeff; fiti.S fito.S fiti.lambdaQ fito.lambdaQ chi];
  catch
  end
end


iOutliers = sum(isoutlier(coeff),2);
avg_coeff = mean(coeff(~iOutliers,:));

[S_i, S_o, lambdaq_i, lambdaq_o, chi] = unpack(avg_coeff);

end











