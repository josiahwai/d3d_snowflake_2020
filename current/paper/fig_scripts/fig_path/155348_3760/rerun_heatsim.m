ccc

shot = 155348;
time_ms = 3760;

lambdaq_i = .0165;  % sol power lengths [m]
lambdaq_x = .006;
lambdaq_o = .003; 
Di = 2;          % Diffusion coeffs [m^2/s]
Dx = 0.1;
Do = 0.2;

% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155348_3760/old_results/eqs.mat')
load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155348_3760/eqs.mat')

sims{1} = heatsim_fit(eqs{1}.gdata, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

plotsim(sims{1})

sims{3} = heatsim_fit(eqs{end}, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

plotsim(sims{3})

save('sims','sims')
