ccc
plotit = 1;
shot = 155354;
time_ms = 3727;
load('/u/jwai/d3d_snowflake_2020/current/dev/3727/eqs.mat')
eq = eqs{end};

root = '/u/jwai/d3d_snowflake_2020/current/';

lambdaq_i = .0063;
lambdaq_o = .006; 
chi_i = .1;
chi_o = .03;

efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];
efit_eq1 = read_eq(shot, time_ms/1000, efit_dir);

load('d3d_obj_mks_struct_6565.mat')

% Load heat flux data q(s,t), s=distance along limiter, and t=time
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
[~,k] = min(abs(t-time_ms));
qperp = qperp(k,:)';

ef = eich_fitter_dev(s', qperp, eq, tok_data_struct, plotit);

% heatsim
lambdaq_i = ef.fiti.lambdaQ / 100
lambdaq_o = ef.fiti.lambdaQ / 100
chi_i = ef.chi_i
chi_o = ef.chi_i
% sim = heatsim_fit(eq, shot, time_ms, lambdaq_i, lambdaq_o, chi_i, chi_o, plotit);


























