addpath(genpath('/home/waij/d3d_snowflake_2020/current'));
ccc

% ========
% SETTINGS
% ========
shot = 155354;
time_ms = 3727;
shotdir = '/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/';

% =======================
% LOAD EQ AND FIT BPROBES
% =======================

efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];
load([shotdir num2str(time_ms) '/eqs.mat'])

eq0 = read_eq(shot, time_ms/1000, efit_dir);
eqf = eqs{end};












% eq_dir = '/fusion/projects/omfit-results/waij/projects/bprobe_errors/EFITtime/OUTPUTS/gEQDSK/';
% eq = read_eq(155354, 4.000, eq_dir);
% struct_to_ws(eq.gdata);
% 
% load('bprobe_data.mat')
% [psi, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, xmp2, ymp2);
% 
% Br = -1./(2*pi*xmp2).*psi_z;
% Bz =  1./(2*pi*xmp2).*psi_r;       
% Bprobe = Br.*cosd(amp2) + Bz.*sind(amp2);
% 
% chisq_expmpi = (fwtmp2 .* (cmpr2 - expmpi)).^2;
% chisq_dB = (fwtmp2 .* (cmpr2 - Bprobe)).^2;
% chisq_Bprobe = (fwtmp2 .* (expmpi - Bprobe)).^2;
% 
% 
% figure
% hold on
% plot(fwtmp2.*expmpi, 'r','linewidth', 1.5)
% plot(fwtmp2.*cmpr2, 'b', 'linewidth', 1.5)
% plot(fwtmp2.*Bprobe, 'g', 'linewidth', 1.5)
% 
% 
% figure
% hold on
% title('Chi^2 Errors')
% p1 = plot(chisq_expmpi,'r');
% p2 = plot(chisq_Bprobe, 'g');
% p3 = plot(chisq_dB, '--k');
% set([p1 p2 p3], 'linewidth', 1.5)
% 
% 
% figure
% hold on
% plot_eq(eq)
% i = 1:length(xmp2);
% scatter(xmp2(i), ymp2(i), 'filled')
% 




























