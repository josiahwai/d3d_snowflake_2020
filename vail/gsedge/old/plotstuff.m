
load eq0.mat
load eq1.mat

eq = eq0;

eq0 = eq;

figure(2)

%.....................
% Compare the profiles

psimag_IN = init.gdata.psimag;
psibry_IN = init.gdata.psibry;

psimag_OT_0 = eq0.psimag;
psibry_OT_0 = eq0.psibry;

psimag_OT_1 = eq2.psimag;
psibry_OT_1 = eq2.psibry;

psi_IN   = linspace(psimag_IN,   psibry_IN,  129);
psi_OT_0 = linspace(psimag_OT_0, psibry_OT_0, 65);
psi_OT_1 = linspace(psimag_OT_1, psibry_OT_1, 65);

psiN_IN   = (psimag_IN - psi_IN)/(psimag_IN - psibry_IN);
psiN_OT_0 = (psimag_OT_0 - psi_OT_0)/(psimag_OT_0 - psibry_OT_0);
psiN_OT_1 = (psimag_OT_1 - psi_OT_1)/(psimag_OT_1 - psibry_OT_1);

% Pressure 

pres_IN   = init.gdata.pres;
pres_OT_0 = eq0.pres;
pres_OT_1 = eq2.pres;

subplot(3,1,1)
plot(psiN_IN, pres_IN, '-b', 'LineWidth', 2)
hold on
plot(psiN_OT_0, pres_OT_0, '--r', 'LineWidth', 2)
hold on
plot(psiN_OT_1, pres_OT_1, '-r', 'LineWidth', 2)

title('Pressure [Pa]')
grid on
xlabel('\psi_N')

% Poloidal current

fpol_IN   = init.gdata.fpol;
fpol_OT_0 = eq0.fpol;
fpol_OT_1 = eq2.fpol;

subplot(3,1,2)
plot(psiN_IN, fpol_IN, '-b', 'LineWidth', 2)
hold on
plot(psiN_OT_0, fpol_OT_0, '--r', 'LineWidth', 2)
hold on
plot(psiN_OT_1, fpol_OT_1, '-r', 'LineWidth', 2)

title('Poloidal Current')
grid on
xlabel('\psi_N')

% Toroidal Current

jtav0 = eq0.jtav;
jtav1 = eq2.jtav;

subplot(3,1,3)
plot(psiN_OT_0, jtav0, '--r', 'LineWidth', 2)
hold on
plot(psiN_OT_1, jtav1, '-r', 'LineWidth', 2)

title('Toroidal Current')
grid on
xlabel('\psi_N')



