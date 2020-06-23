% Sometimes efit mis-predicts the snowflake type vs IR.
% Scan through g-files and IR to find an example shot/time where this
% occurs. 

close all;
shot = 155354;

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));
warning('off','all')

eq_dir = '/u/jwai/d3d_snowflake_2020/current/inputs/eqs/cake/155354/';

% efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];

% find available times
d = dir(eq_dir);
times = [];
for k = 3:length(d)
  times = [times; str2num(d(k).name(end-3:end))];
end

% analyze the snowflake types
for iTime = 1:length(times)
%   try
  time_ms = times(iTime);
  
  % ------- 
  % EFIT eq
  % -------
  eq = read_eq(shot, time_ms/1000, eq_dir);
  snow = analyzeSnowflake(eq.gdata); 
  
  % -------  
  % IRTV 
  % -------
  
  % Load heat flux data q(s,t), s=distance along limiter, and t=time
  qperp_dir  = [root 'inputs/qperp/'];
  qperp_data = ['qperp_' num2str(shot) '.mat'];
  
  load([qperp_dir qperp_data])  % loads q, s, and t
  
  [~,k] = min(abs(t-time_ms));
  qir = qperp(k,:)'/100;
  s = s/100;
  
  pkthresh = 2*median(qir);
  
  iX = find(s>1.15 & s<1.45);
   
  [qirmaxX] = findpeaks(qir(iX),'NPeaks',1,'sortstr','descend',...
    'minpeakheight', pkthresh, 'minpeakprominence', pkthresh);
  
  snowPlus = false;
  if isempty(qirmaxX), snowPlus=true; end
  
  % mismatch?
%   if snow.snowPlus ~= snowPlus
%     disp([shot time_ms])
    
    plot_eq(eq,18)
    title(num2str(time_ms))
    set(gcf,'position', [800 98 329 312])
    
    
    figure(201)
    plot(s,qir)  
    title(num2str(time_ms))
    set(gcf,'position', [776 481 443 219])
    
    
    clf(18)
    clf(201)
    
%   end  
   
%   catch
%   end
end
































