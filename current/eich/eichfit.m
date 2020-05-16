shot = 155354;
root = '/u/jwai/d3d_snowflake_2020/current/';

% load heat flux
qperp_dir  = [root 'inputs/qperp/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];
load([qperp_dir qperp_data])  % loads q, s, and t
[~,iTime] = min(abs(t - median(t_3pks)));
time_ms = t(iTime);
q = qperp(iTime,:)';

% load eq and calculate flux expansion
eqdir = [root 'inputs/eqs/cake/' num2str(shot)];
eq = read_eq(shot, time_ms, eqdir);
snow = analyzeSnowflake(eq);
[fexpi,fexpo] = unpack(calc_fluxexp(eq.gdata, snow.rSPP, snow.zSPP));


% remove gap from limiter
s = s';
iGap = find(s < 170,1,'last');
dgap = s(iGap+1) - s(iGap);
s(iGap(end)+1:end) = s(iGap(end)+1:end) - dgap;
ii = find(s<120);
io = find(s>145);


% Fit eich profile to the outer peak
eichfun_o = @(q0, S, lambdaQ, fExp, s0, qBg, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBg;

ft = fittype(eichfun_o);
options = fitoptions(ft);
options.StartPoint = [75 1 1 fexpo 150 4];
options.Lower = [0 0 0 fexpo 0 0];  % lock flux expansion
options.Upper = [inf inf fexpo inf inf];


fito = fit(s(io), q(io), ft, options); 
q0      = fito.q0;
S       = fito.S;
lambdaQ = fito.lambdaQ;
fExp    = fito.fExp;
s0      = fito.s0;
qBg     = fito.qBg;

sfito = linspace(0,200,1000);
qfito = eichfun_o(q0,S,lambdaQ,fExp,s0,qBg,sfito);
qfito(qfito<1.03*qBg) = nan;


% Fit eich profile to the inner peak
eichfun_i = @(q0, S, lambdaQ, fExp, s0, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (-(x - s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x - s0))./S);

ft = fittype(eichfun_i);
options = fitoptions(ft);
options.StartPoint = [75 1 1 fexpi 100]; 
options.Lower = [0 0 0 fexpi 0];
options.Upper = [inf inf fexpi inf];

fiti = fit(s(ii), q(ii), ft, options); 
q0      = fiti.q0;
S       = fiti.S;
lambdaQ = fiti.lambdaQ;
fExp    = fiti.fExp;
s0      = fiti.s0;

sfiti = linspace(0,200,1000);

qfiti = eichfun_i(q0,S,lambdaQ,fExp,s0,sfiti);
qfiti(qfiti<1.03*qBg) = nan;


% plot the fits
figure
hold on
plot(s, q, 'k')
plot(sfiti, qfiti, 'r', 'linewidth', 1.5)
plot(sfito, qfito, 'r', 'linewidth', 1.5)


% Calculate thermal diffusivity from the outer strike point fit
Ti = 50 * 11600;  % [K]
Te = 50 * 11600;  % [K]
k = 1.38e-23;     % Boltzmann constant [(m^2*kg)/(s^2*K)]
mi = 3.344e-27;   % Deuterium mass [kg]
cs = sqrt(k*(Ti + Te)/mi); % sound speed
L = 5;                     % connection length varies from 1.5 to ~20
tau = L/cs;                % time of flight
chiTemp = (fito.S/100)^2/tau;   % thermal diffusivity in units [m^2/s]


% Compute the flux gradient at strike point
struct_to_ws(eq.gdata);
[~, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, snow.rSPP(2), snow.zSPP(2));

gradpsi = sqrt(psi_r^2 + psi_z^2);
chi = chiTemp*gradpsi^2; % thermal diffusivity in units [Wb^2/s]


fprintf("chi = %3.3f Wb^2/s\n", chi)
fprintf("lambdaQ = %3.5f m\n", fito.lambdaQ/100)

















