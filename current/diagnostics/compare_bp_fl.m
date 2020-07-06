addpath(genpath('/home/waij/d3d_snowflake_2020/current'));
ccc

% ========
% SETTINGS
% ========
shot = 155330;
time_ms = 3660;
shotdir = ['/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/' num2str(shot) '_sfm/'];

% =======================
% LOAD EQ AND FIT BPROBES
% =======================

efit_dir = ['/p/omfit/users/jwai/projects/snowflake_efits/EFITtime/OUTPUTS/' ...
  num2str(shot) '/gEQDSK/'];
load([shotdir num2str(time_ms) '/eqs.mat'])

% eq0 = read_eq(shot, time_ms/1000, efit_dir);
eq0 = eqs{2};
eqf = eqs{end};

load(['bp_fl_data' num2str(shot)])
bp_fl = getfield(bp_fl_data, ['t' num2str(time_ms)]);
struct_to_ws(bp_fl);
struct_to_ws(bp_fl_data.template);

fwtmp2 = fwtmp2 / max(fwtmp2);
% fwtmp2 = FWTMP2;


Bp0 = calc_mag_probe(eq0, XMP2, YMP2, AMP2);
Bpf = calc_mag_probe(eqf, XMP2, YMP2, AMP2);

chisq0 = (fwtmp2 .* (expmpi - Bp0)).^2;  % fitting w
chisqf = (fwtmp2 .* (expmpi - Bpf)).^2;


% ========
% PLOT IT
% ========
figure
subplot(2,1,1)
hold on
title('FitWeight x Bp')
plot( fwtmp2.*expmpi, '--k')
plot( fwtmp2.*Bp0, 'r')
plot( fwtmp2.*Bpf, 'b')
xline(48,'--k');
xline(51,'--k');
legend('Measured', 'Calculated: EFIT', 'Calculated: EFIT+IRTV','location','southwest')


subplot(2,1,2)
hold on
title('Chi^2 B-Probes')
plot( chisq0, 'r')
plot( chisqf, 'b')
xline(48,'--k');
xline(51,'--k');
legend('EFIT', 'EFIT+IRTV')

set(gcf,'position',[440 177 657 526])


plot_eq(eq0)
hold on
i = 48:51;
scatter(XMP2(i), YMP2(i), 'filled')
set(gcf,'position',[740 147 418 553])






































