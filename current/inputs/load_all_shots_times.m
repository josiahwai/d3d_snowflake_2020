% CREATE LIST OF SHOTS/TIMES WHERE DATA EXISTS FROM: EFIT, CAKE, AND IRTV

shotlist = 155328:155355;

root = '/u/jwai/d3d_snowflake_2020/current/';
all_shots_times = [];

for shot = shotlist
  
  % times with cake
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  d = dir(cake_dir);
  t_cake = [];
  for k = 3:length(d)
    t_cake = [t_cake; str2num(d(k).name(end-3:end))];
  end
  
  % times with efit
  efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
    num2str(shot) '/gEQDSK/'];
  
  d = dir(efit_dir);
  t_efit = [];
  for k = 3:length(d)
    t_efit = [t_efit; str2num(d(k).name(end-3:end))];
  end
  
  
  % ir times
  qperp_dir  = [root 'inputs/qperp/'];
  qperp_data = ['qperp_' num2str(shot) '.mat'];
  load([qperp_dir qperp_data])  % loads q, s, and t
  
  t_ir = [t_2pks t_3pks];
  
  % times with efit and cake
  t_eq = intersect(t_cake, t_efit);
  
  % times with IRTV within 10ms of efit/cake times
  [dt,k] = min(abs(t_eq - t_ir));
  t_good = t_eq(k(dt<10));  % within 10 ms
  
  shots = shot * ones(size(t_good));
  
  all_shots_times = [all_shots_times; shots t_good];

end

fn = '/u/jwai/d3d_snowflake_2020/current/inputs/all_shots_times';

save(fn, 'all_shots_times')









