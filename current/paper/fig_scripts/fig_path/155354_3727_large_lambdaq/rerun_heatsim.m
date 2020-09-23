ccc

load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155354_3727_large_lambdaq/eqs.mat')

shot = 155354;
time_ms = 3727;

lambdaq_i = .009;  % sol power lengths [m]
lambdaq_x = .006;
lambdaq_o = .002; 
Di = 2;          % Diffusion coeffs [m^2/s]
Dx = 0.1;
Do = 0.15;

sims{1} = heatsim_fit(eqs{1}.gdata, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

sims{1}.qO = smooth(sims{1}.qO, 3);

plotsim(sims{1})


lambdaq_i = .009;  % sol power lengths [m]
lambdaq_x = .006;
lambdaq_o = .0024; 
Di = 2;            % Diffusion coeffs [m^2/s]
Dx = 0.1;
Do = 1.2;

sims{3} = heatsim_fit(eqs{end}, shot, time_ms, lambdaq_i, lambdaq_x, lambdaq_o, ...
  Di, Dx, Do);

plotsim(sims{3})

save('sims', 'sims')
















