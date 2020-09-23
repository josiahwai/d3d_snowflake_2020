ccc

load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155350_3893/eqs.mat')

shot = 155350;
time_ms = 3893;

lambdaq_i = .013;  % sol power lengths [m]
lambdaq_x = .003;
lambdaq_o = .003; 
Di = 1.8;          % Diffusion coeffs [m^2/s]
Dx = 0.05;
Do = 0.4;

sims{1} = heatsim_fit(eqs{1}.gdata, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

sims{1}.qO = smooth(sims{1}.qO, 3);

plotsim(sims{1})


lambdaq_i = .013;  % sol power lengths [m]
lambdaq_x = .003;
lambdaq_o = .003; 
Di = 1.8;            % Diffusion coeffs [m^2/s]
Dx = 0.05;
Do = 0.4;

sims{3} = heatsim_fit(eqs{end}, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

plotsim(sims{3})

save('sims', 'sims')
















