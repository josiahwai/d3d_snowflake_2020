load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfm.mat')
sim1 = sim;
load('/u/jwai/d3d_snowflake_2020/current/paper/analysis/load_sims_data/sims_sfp_sp.mat')


djmax = [sim.djmax sim1.djmax];
j0 = [sim.j0max sim1.j0max];


mean(djmax./j0)


mean(abs(djmax./j0))


