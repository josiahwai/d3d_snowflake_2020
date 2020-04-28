% Set up shared variables with OUTFUN
clear spec init config gsdesign
history.x = [];
history.fval = [];

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[iJob,shot,time_ms,rxP,rxS,zxP,zxS] = unpack(dum.args);
x0 = [rxP rxS zxP zxS];

fprintf('\n\nshot: %d time: %d \n\n', shot, time_ms);

lb = x0-0.12;
ub = x0+0.12;

options = optimoptions(@fmincon, ... %'OutputFcn',@outfun,...
    'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',...
    0.002,'DiffMinChange',0.002);

xsol = fmincon(@(xp)sfd_costfunction(xp,shot,time_ms),x0,[],[],[],[],lb,ub,[],options);








      
      











