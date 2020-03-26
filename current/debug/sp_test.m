shot = 155355; 
time_ms = 4400;
root = '/u/jwai/d3d_snowflake_2020/current/';

opts.plotit = 1;
opts.saveit = 1;
opts.root = root;
opts.iSim = 1;

efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
efit_eq = read_eq(shot, time_ms/1000, efit_dir);
opts.saveDir = '/u/jwai/d3d_snowflake_2020/current/debug/';
if ~exist(opts.saveDir, 'dir')
  mkdir(opts.saveDir);
end
heatsim_batch2(efit_eq.gdata, shot, time_ms, opts);

% J = heatsim_cost(efit_eq.gdata,shot,time_ms,opts);