% -------
% TESTING
% -------
% clear; clc; close all;
% shot = 165288;
% time_ms = 1500;
% root = '/u/jwai/d3d_snowflake_2020/current/';
% efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
% eq = read_eq(shot, time_ms/1000, efit_dir);
% x = featureMapSnow(eq,shot,time_ms)

% ----------------------
% FEATURE MAP SNOW
% ----------------------
function x = featureMapSnow(eq,shot,time_ms)

plotit = 0;

if isfield(eq,'gdata'), eq = eq.gdata; end

x1 = analyzeSnowflake(eq,plotit);
x2 = heatsim_ml(eq,shot,time_ms,plotit);

struct_to_ws(x1);
struct_to_ws(x2);

x = [rx zx psix sSP_heat s_qmax s_qirmax qmaxN ...
  qirmaxN Apk_qmaxN Apk_qirmaxN psi_qmax psi_qirmax, npeaks_q, npeaks_qir];

x(isnan(x)) = 0;














