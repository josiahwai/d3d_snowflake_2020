clear all; clc; close all; warning('off','all');

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);

% =========================
% LOAD AND SIMULATE EFIT EQS
% =========================

load('eq0')
N = 16;
th = linspace(0,2*pi,N+1);
th(end) = [];

snow0 = analyzeSnowflake(eq0);
xp0 = [snow0.rx snow0.zx];
sp0 = [snow0.rSPP snow0.zSPP];

sps = []; 
xps = [];
eqs = {};
for k = 1:length(th)
  
%   dxp(k,:) =  0.02 * [0 cos(th(k)) 0 sin(th(k))];
  dxp(k,:) = k / length(th) * .04 * [0 0 0 1];
  xps(k,:) = xp0 + dxp(k,:);
  
  eqs{k} = designeq_ml( xps(k,:), shot, time_ms);
  
  snow = analyzeSnowflake(eqs{k});
  
  sps(k,:) = [snow.rSPP snow.zSPP];
  dsp(k,:) = sps(k,:) - sp0;
end



% plot it 
figure
c = flip(cool);
N = length(th);

for k = 1:N
  
  co = c(floor(k/N*length(c)),:);
  
  subplot(1,2,1)
  hold on
  scatter( xps(k,1:2), xps(k,3:4), 40, 'filled', 'markerfacecolor', co)
  axis([1.0456    1.1790   -1.4741   -1.2741])
  
  subplot(1,2,2)
  hold on
  scatter( sps(k,1:2), sps(k,3:4), 40, 'filled', 'markerfacecolor', co)
  scatter( sp0(1:2), sp0(3:4), 40, 'k', 'filled')
  axis([1.2782    1.2988   -1.3636   -1.3628])
  pause
end


% plot limiter and strike points
subplot(1,2,2)
load('d3d_obj_mks_struct_6565.mat')
rlim = tok_data_struct.limdata(2,:); 
zlim = tok_data_struct.limdata(1,:);
plot(rlim, zlim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
axis([1 1.4 -1.4 -1])
axis([1.2782    1.2988   -1.3636   -1.3628])




















