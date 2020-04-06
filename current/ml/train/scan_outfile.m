shot = 155328;
time_ms = 3850;
ic = 1;

clc
root = '/u/jwai/d3d_snowflake_2020/current/';
jobdir = [root 'ml/train/jobs/' num2str(shot) '/' num2str(time_ms) '/' num2str(ic) '/'];
d = dir([jobdir '*.out']);
fn = d.name;
fid = fopen([jobdir fn]);
tline = fgetl(fid);

cost = [];
xp = [];
k = 0;
while ischar(tline)
  if contains(tline,'Cost:')    
    c = textscan(tline, '%s %f %s %f %f %f %f'); 
    k = k+1;
    [~,cost(k),~,rxP,rxS,zxP,zxS] = unpack(c);
    xp(k,:) = [rxP rxS zxP zxS];
    
%     if ~contains(fgetl(fid), 'read_eq')
%       k = k+1;
%       [~,cost(k),~,rxP,rxS,zxP,zxS] = unpack(c);
%       xp(k,:) = [rxP rxS zxP zxS];
%     end
  end  
  tline = fgetl(fid);
end

% % close all
% efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
% efit_eq = read_eq(shot, time_ms/1000, efit_dir);
% plot_eq(efit_eq)
% hold on
% axis equal
% axis([1.0 1.5 -1.4 -0.9])
plot(xp(:,1:2), xp(:,3:4),'b','linewidth',1.5)
scatter(xp(end,1:2), xp(end,3:4),'r','filled')

[cost' xp]



















