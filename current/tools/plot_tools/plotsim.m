function sim = plotsim(sim)

figure(27)
hold on
struct_to_ws(sim);

% normalize heat flux
load('d3d_obj_mks_struct_6565.mat')
limdata = tok_data_struct.limdata;
slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
[rir,zir] = calcLimDistanceInv( slimtot - sir, limdata);

smax = 1.61;
qir(sir>smax) = 0;

iI = find(sir < 1.2);
iX = find(sir > 1.2 & sir < 1.4);
iO = find(sir > 1.4);

Pmeas1 = trapz(sir(iI), 2*pi*rir(iI).*qir(iI)');
Pmeas2 = trapz(sir(iX), 2*pi*rir(iX).*qir(iX)');
Pmeas4 = trapz(sir(iO), 2*pi*rir(iO).*qir(iO)');
Pmeas = Pmeas1 + Pmeas2 + Pmeas4;

Psim1 = abs(trapz(sI, 2*pi*r_qmax(1)*qI));
if ~isnan(r_qmax(2))
  Psim2 = abs(trapz(sX, 2*pi*r_qmax(2)*qX));
else
  Psim2 = 0;
end
qO_dum = qO;
qO_dum(sO>smax) = 0;
Psim4 = abs(trapz(sO, 2*pi*r_qmax(3)*qO_dum));
Psim = Psim1 + Psim2 + Psim4;

qirN = qir / Pmeas;

% qIN = qI * Pmeas1 / Psim1 / Pmeas;
qIN = qI * max(qir(iI)) / max(qI) / Pmeas;
qXN = qX * (Pmeas2 + Pmeas4) / (Psim2 + Psim4) / Pmeas;
qON = qO * (Pmeas2 + Pmeas4) / (Psim2 + Psim4) / Pmeas;

qN = [qIN; qXN; qON];


% plot normalized heat flux
plot(sir, qirN, 'k', 'linewidth', 1.5)

s = [sI sX sO];
plot(s, qN, 'linewidth', 1.5)

axis([0.8 1.8 0 1.05* max( max(qirN), nanmax(qN))])
set(gcf,'position', [718 -93 535 190])


sim.qirN = qirN;
sim.qIN = qIN;
sim.qON = qON;
sim.qXN = qXN;
sim.qN = qN;
sim.s = s;


















