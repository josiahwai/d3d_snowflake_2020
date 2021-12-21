% topdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/cakes_early_t/';
topdir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/cakes_early_t/';


d = dir(topdir);
d(1:2) = [];

shots = [];
times = [];
shots(1) = str2num( d(1).name(2:7));
times(1) = str2num( d(1).name(end-4:end));

for i = 2:length(d)
  shots(i) = str2num( d(i).name(2:7));
  times(i) = str2num( d(i).name(end-4:end));
  
  if shots(i) ~= shots(i-1)
    clf
    shots(i)
    times(i)
    eq = read_eq( shots(i), times(i), topdir);
    plot_eq(eq);
    set(gcf,'position',[735 302 388 298])
    title( [ num2str(shots(i)) ': ' num2str(times(i)) 'ms'])
  end
end


save('shots','shots');
save('times','times');



  
  
  
  