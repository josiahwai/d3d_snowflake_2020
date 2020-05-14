close all;
shot = 155354;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));
warning('off','all')

% efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];

efit_dir = '/p/omfit/users/jwai/projects/155354_efits/EFITtime/OUTPUTS/gEQDSK';

% find times
d = dir(efit_dir);
times = [];
for k = 3:length(d)
  times = [times; str2num(d(k).name(end-3:end))];
end


for iTime = 1:length(times)
  
  time_ms = times(iTime);
  efit_eq = read_eq(shot, [time_ms+5 time_ms-5]/1000, efit_dir);
  clf
  plot_eq(efit_eq);
  axis([1 1.5 -1.4 -0.9])
  title(num2str(time_ms))
  pause
end
  
  
  
  
  
  
  
  
  
  
  
  