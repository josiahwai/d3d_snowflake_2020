d = dir(pwd);
d(1:2) = [];

for k = 1:length(d)
  
  t = str2num(d(k).name(end-3:end));
  
  if ~ismember(t,t_sim)
    movefile(d(k).name, './old')
  end
end










