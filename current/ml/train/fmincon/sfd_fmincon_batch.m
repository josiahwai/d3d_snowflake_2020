% Set up shared variables with OUTFUN
clear spec init config gsdesign
history.x = [];
history.fval = [];

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

dum = load('args');
[shot,time_ms,rxP,rxS,zxP,zxS] = unpack(dum.args);
xp0 = [rxP rxS zxP zxS];

fprintf('\n\nshot: %d time: %d \n\n', shot, time_ms);

lb = xp0-0.12;
ub = xp0+0.12;

options = optimoptions(@fmincon, ... %'OutputFcn',@outfun,...
    'Display','iter','Algorithm','interior-point','FiniteDifferenceStepSize',...
    0.005,'DiffMinChange',0.005);

xsol = fmincon(@(xp)sfd_costfunction(xp,shot,time_ms),xp0,[],[],[],[],lb,ub,[],options);

% run a final time and save it
sfd_costfunction(xsol,shot,time_ms,1)





      
      











