
shotlist = [155328:155355 155478 165288];
savedir = '/u/jwai/d3d_snowflake_2020/current/inputs/times/';

for i = 1:length(shotlist)
  shot = shotlist(i);
  
  % times with cake
  cake_dir = [root 'inputs/eqs/cake/' num2str(shot)];
  d = dir(cake_dir);
  t = [];
  for k = 3:length(d)
    t = [t; str2num(d(k).name(end-3:end))];
  end
  
  fn = [savedir 'times_' num2str(shot) '.mat'];
  save(fn, 't')
end
