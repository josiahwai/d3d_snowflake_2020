
ccc
shotlist = [ 155328 155329 155330 155331 155332 155340 ];


P2_all = [];
P4_all = [];
drsplit_all = [];
lambdaqo_all = [];
lambdaqi_all = [];

for ishot = 1:length(shotlist)
  shot = shotlist(ishot);
  load(['sim' num2str(shot) '.mat'])
  struct_to_ws(sim);
  
  P2_all = [P2_all P2];
  P4_all = [P4_all P4];
  drsplit_all = [drsplit_all drsplit];
  lambdaqo_all = [lambdaqo_all lambdaq_o];
  lambdaqi_all = [lambdaqi_all lambdaq_i];
end

clear sim
sim.P2 = P2_all;
sim.P4 = P4_all;
sim.drsplit = drsplit_all;
sim.lambdaq_o = lambdaqo_all;
sim.lambdaq_i = lambdaqi_all;

savedir = '/u/jwai/d3d_snowflake_2020/current/lambdaq_split/';
save([savedir 'sim_all'], 'sim');
  
fit_lambda_q
  
