dum = load('args');
[shot,time_ms,constrain_sp] = unpack(dum.args);


load(['bp_fl_data' num2str(shot)])
bp_fl = getfield(bp_fl_data, ['t' num2str(time_ms)]);
struct_to_ws(bp_fl);
struct_to_ws(bp_fl_data.template);

fwtmp2 = fwtmp2 / max(fwtmp2);
% fwtmp2 = FWTMP2;

load('eq1.mat')
load('eqf.mat')
Bp0 = calc_mag_probe(eq1, XMP2, YMP2, AMP2);
Bpf = calc_mag_probe(eqf, XMP2, YMP2, AMP2);

chisq0 = (fwtmp2 .* (expmpi - Bp0)).^2;  % fitting w
chisqf = (fwtmp2 .* (expmpi - Bpf)).^2;


% ========
% PLOT IT
% ========
figure
subplot(2,1,1)
hold on
title('FitWeight * B')
plot( fwtmp2.*expmpi, '--k', 'linewidth', 1)
plot( fwtmp2.*Bp0, 'r', 'linewidth', 1)
plot( fwtmp2.*Bpf, 'b', 'linewidth', 1)
ylabel('Field [T]')
xlabel('Bprobe Index')
% xlim([42 62])
legend('Measured', 'EFIT', 'EFIT+IRTV','location','southwest')


subplot(2,1,2)
hold on
hold on
plot( fwtmp2.*(expmpi - Bp0), 'r', 'linewidth', 1)
plot( fwtmp2.*(expmpi - Bpf), 'b', 'linewidth', 1)
title('FitWeight * (Calculated - Measured)')
ylabel('Field [T]')
xlabel('Probe index')


% title('Chi^2 B-Probes')
% plot( chisq0, 'r', 'linewidth', 1.5)
% plot( chisqf, 'b', 'linewidth', 1.5)
% ylabel('Chi^2')
% xlabel('Bprobe Index')
% load('chisq_smooth.mat')
% plot( sum(chisq_smooth), '--k', 'linewidth', 1.5)
% xlim([42 62])
legend('EFIT', 'EFIT+IRTV')

set(gcf,'position',[440 177 657 526])


% plot_eq(eq0)
% hold on
% i = 48:51;
% scatter(XMP2, YMP2, 'b', 'filled')
% scatter(XMP2(i), YMP2(i), 'r', 'filled')
% set(gcf,'position',[740 147 418 553])







% figure()
% hold on
% plot( fwtmp2.*(expmpi - Bp2), 'r', 'linewidth', 1.5)
% plot( fwtmp2.*(expmpi - Bpf), 'b', 'linewidth', 1.5)
% title('Fit weight * (Calculated - Measured)')
% ylabel('B [T]')
% xlabel('Probe index')
% set(gcf,'position',[437 406 635 294])
% plot( fwtmp2.*(Bp2 - Bpf), 'b', 'linewidth', 1.5)































