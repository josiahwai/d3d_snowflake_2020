shot = 155354;

% load heat flux
root = '/u/jwai/d3d_snowflake_2020/current/';
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t

% remove gap from limiter
s = s';
iGap = find(s < 170,1,'last');
dgap = s(iGap+1) - s(iGap);

s_nogap = s;
s_nogap(iGap(end)+1:end) = s_nogap(iGap(end)+1:end) - dgap;

load('d3d_obj_mks_struct_6565.mat');
limdata = tok_data_struct.limdata;

slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);


[r,z] = calcLimDistanceInv(slimtot - s/100, limdata);

s1 = slimtot - calcLimDistance(r,z, limdata);


sgap = 1.726:0.01:1.726+dgap/100;
[rgap,zgap] = calcLimDistanceInv(slimtot - sgap, limdata);

figure
hold on
plot(r,z,'r')
plot(rgap,zgap,'b','linewidth', 2)
























