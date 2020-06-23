% Fits the eich profile to the heat flux data and finds the strike points
%         Takes the shadowing effect into account by cutting off a larger
%         portion of s
%
% INPUTS: s - distance along limiter, qperp - heat flux profile versus 
%         position (s), eq - equilibrium, tok_data_struct - tokamak
%         geometry struct, plotit - flag to plot some results
%
% OUTPUTS: rsp, zsp - coordinates of strike points as predicted by the heat
%         flux, ssp - position along limiter of strike pts, fiti / fito -
%         inner and outer fits for the eich profile with coefficient
%         information as defined in ref. 
%
% ref: T. Eich, "Scaling of the tokamak near the scrape-off layer H-mode 
%      power width and implications for ITER"


function eichfit = eich_fitter_geo(s, qperp, eq, ...
  tok_data_struct, plotit)

if ~exist('plotit','var'), plotit = 0; end
if isfield(eq,'gdata'), eq = eq.gdata; end

% define the eich functions
eichfun_i = @(q0, S, lambdaQ, fExp, s0, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (-(x - s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x - s0))./(S*fExp));

eichfun_o = @(q0, S, lambdaQ, fExp, s0, qBg, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./(S*fExp)) + qBg;


% determine flux expansion from geometry
snow = analyzeSnowflake(eq);
[fexpi,fexpo] = unpack(calc_fluxexp(eq, snow.rSPP([1 end]), snow.zSPP([1 end])));
    

% index of heat flux regions (inner, outer, middle)
if size(s,2) ~= 1, s = s'; end

s_shelftop = 184.7;
s_shelfbot = 163.8;
ds = s_shelftop - s_shelfbot;

iuse = s <= s_shelfbot | s >= s_shelftop;

s_geo = s;
q_geo = qperp;
s_geo(~iuse) = [];
for i = 1:length(s_geo)
  if s_geo(i) >= s_shelftop
    s_geo(i) = s_geo(i) - ds;
  end
end
q_geo(~iuse) = [];

ii = find(s_geo < 120);
io = find(s_geo > 145);
ix = find( s_geo >= 120 & s_geo <= 145);


% fit the outer peak to find strike pt
% ....................................

ft = fittype(eichfun_o);
options = fitoptions(ft);

options.StartPoint = [max(q_geo(io)) 0.2 0.3 fexpo 160 4];
options.Lower = [0 0 0 fexpo 0 0]; 
options.Upper = [inf inf inf fexpo inf inf];  % lock flux expansion

[fito, gofo] = fit(s_geo(io), q_geo(io), ft, options);
sspo =  fito.s0;


% fit the middle peak, if the peak is there
% .........................................
pkthresh = 3*median(qperp);

[qpkx,ipkx] = findpeaks(q_geo(ix), 'NPeaks',1,'sortstr','descend',...
  'minpeakheight', pkthresh, 'minpeakprominence', pkthresh);

% estimate strike point as location where qperp = 3/4*q_peak on inbd side
[~,k] = min(abs(q_geo(ix(1:ipkx)) - 0.85*qpkx));
sspx = s(ix(k));

if isempty(sspx), sspx = nan; end


% fit the inner peak to find strike pt
% ....................................
ft = fittype(eichfun_i);
options = fitoptions(ft);
options.StartPoint = [max(q_geo(ii)) 1 0.6 fexpi 100];
options.Lower = [0 0 0 fexpi 0];
options.Upper = [inf inf inf fexpi inf];

[fiti, gofi] = fit(s(ii), qperp(ii), ft, options);
sspi =  fiti.s0;

% convert strike points from s to (r,z)
ssp = [sspi sspx sspo]/100;

limdata = tok_data_struct.limdata;
slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
[rsp,zsp] = calcLimDistanceInv(slimtot - ssp, limdata);


% Calculate thermal diffusivity (related to S param in Eich model)
% ................................................................
Ti = 50 * 11600;  % [K]
Te = 50 * 11600;  % [K]
k = 1.38e-23;     % Boltzmann constant [(m^2*kg)/(s^2*K)]
mi = 3.344e-27;   % Deuterium mass [kg]
cs = sqrt(k*(Ti + Te)/mi); % sound speed

% Outer chi
L = 3;                          % connection length [m]
tau = L/cs;                     % time of flight
chiTemp = (fito.S/100)^2/tau;   % thermal diffusivity in units [m^2/s]    
struct_to_ws(eq);               % Compute the flux gradient at strike point
[~, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, snow.rSPP(end), snow.zSPP(end));
gradpsi = sqrt(psi_r^2 + psi_z^2);
chi_o = chiTemp*gradpsi^2;        % thermal diffusivity in units [Wb^2/s]

% Inner chi
L = 2;                          % connection length [m]
tau = L/cs;                     % time of flight
chiTemp = (fiti.S/100)^2/tau;   % thermal diffusivity in units [m^2/s]    
struct_to_ws(eq);               % Compute the flux gradient at strike point
[~, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, snow.rSPP(1), snow.zSPP(1));
gradpsi = sqrt(psi_r^2 + psi_z^2);
chi_i = chiTemp*gradpsi^2;        % thermal diffusivity in units [Wb^2/s]


% write qperp peak info
if ~isnan(sspx)
  qpkx = max(qperp(ix));
else
  qpkx = nan;
end
qpks = [max(qperp(ii)) qpkx max(qperp(io))];


eichfit = struct('rsp', rsp, 'zsp', zsp, 'ssp', ssp, 'chi_i', chi_i, ...
  'chi_o', chi_o, 'qpks', qpks, 'fiti', fiti, 'fito', fito, 'gofi', ...
  gofi, 'gofo', gofo);


if plotit
  figure(199)
  clf
  hold on
  plot( s, qperp,'.k')
  
  % plot inner
  plot( s_geo(ii), fiti(s_geo(ii)), 'r', 'linewidth', 1.5)
 
  % restore the gap and plot outer
  q_fito = fito(s_geo(io));
  s_fito = s_geo(io);  
  for i = 1:length(s_fito)
    if s_fito(i) > s_shelfbot
      s_fito(i) =s_fito(i) + ds;
    end
  end
  
  plot( s_fito, q_fito, 'g', 'linewidth', 1.5)   
  
  xline(sspi);
  if ~isnan(sspx), xline(sspx, 'b', 'linewidth', 2); end
  xline(sspo);
  
%   add_lim_colorcode(199)
%   xlim([84 213])
%   set(gcf,'position', [680 608 560 272])
end

end









