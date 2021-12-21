clear all; clc; close all; warning('off','all');

wts = logspace(1,4.7,100);

root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));

shot = 155350;
time_ms = 3893;
constrain_sp = 0;

% Weight scan
load('xps155350_3893');
xpf = xps{end};

% ==========================
% LOAD AND SIMULATE EFIT EQS
% ==========================

% load and simulate
efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];

% efit01
efit_eq0 = read_eq(shot, time_ms/1000, efit_dir);
efit_snow0 = analyzeSnowflake(efit_eq0);
efit_xp0 = [efit_snow0.rx efit_snow0.zx];

% hybrid eq
clear iter init spec config gsdesign 
efit_eq1 = designeq_bp(efit_xp0,shot,time_ms,[],1e5);
efit_eq1 = rmfield(efit_eq1, {'r', 'b', 'p'});

% B-probe fit
load(['bp_fl_data' num2str(shot)])
bp_fl = getfield(bp_fl_data, ['t' num2str(time_ms)]);
struct_to_ws(bp_fl);
struct_to_ws(bp_fl_data.template);

fwtmp2 = fwtmp2 / max(fwtmp2);
Bp = calc_mag_probe(efit_eq1, XMP2, YMP2, AMP2); 


% initialize
eqs = {efit_eq1};
chisq = {(fwtmp2 .* (expmpi - Bp)).^2};
ssq = {sum(chisq{1})};
Bps = {Bp};
snows = {analyzeSnowflake(eqs{1})};





for k = 2:length(wts)+1
  try
    eq = designeq_bp(xpf,shot,time_ms,[],wts(k+1));
    
    eqs{k} = rmfield(eq, {'r', 'b', 'p'});
    Bps{k} = calc_mag_probe(eqs{k}, XMP2, YMP2, AMP2);
    chisq{k} = (fwtmp2 .* (expmpi - Bps{k})).^2;
    ssq{k} = sum(chisq{k});
    snows{k} = analyzeSnowflake(eqs{k});
  catch
    warning('on','all')
    warning(['Didnt run at ' num2str(k)]);
  end
end

bprobe_weight_scan_data = struct('eqs',eqs,'Bps',Bps,'chisq',chisq,...
  'ssq',ssq,'snows',snows, 'wts', wts);


fn = ['/u/jwai/d3d_snowflake_2020/current/diagnostics/bprobe_weight_scan_data'...
  num2str(shot) '.mat'];
save(fn, 'bprobe_weight_scan_data')






























