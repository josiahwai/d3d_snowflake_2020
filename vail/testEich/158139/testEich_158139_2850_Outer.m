
% Load tokamak definition

load('/u/pvail/d3d_snowflake_2019/invMap/d3d_obj_mks_struct_129129.mat')

limdata = tok_data_struct.limdata;

% Configure the plots

plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([1.0 1.5 -1.4 -0.9])
xlabel('R [m]')
ylabel('Z [m]')
title([int2str(158139) ': ' int2str(2850) ' ms'])
    
% Load the equilibrium

eq = read_eq(158139, 2.85, '/u/pvail/d3d_snowflake_2019/testEich/158139');

rg = eq.gdata.rg; 
zg = eq.gdata.zg;

zmaxis = eq.gdata.zmaxis;

bzero = eq.gdata.bzero;
rzero = eq.gdata.rzero;

psizr  = eq.gdata.psizr;
psibry = eq.gdata.psibry;

[psizr257257, rg257, zg257] = regrid(rg, zg, psizr, 257, 257);

% Compute Br and Bz at the grid points

[dpsidr, dpsidz] = rzGrad(psizr257257, rg257, zg257);

Br257257 = zeros(size(psizr257257,1),size(psizr257257,2));
Bz257257 = zeros(size(psizr257257,1),size(psizr257257,2));

for ii = 1:length(rg257)
    for jj = 1:length(zg257)
        
        Br257257(jj,ii) = -1/(2*pi*rg257(ii))*dpsidz(jj,ii);
        Bz257257(jj,ii) =  1/(2*pi*rg257(ii))*dpsidr(jj,ii);
        
    end
end

% Analyze the equilibrium

rxPL =  1.20;
zxPL = -1.20;

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr257257, rxPL, zxPL, rg257, zg257);

contour(rg257, zg257, psizr257257, [psixPL psixPL], 'k', 'LineWidth', 2)

plot(rxPL, zxPL, 'xk', 'Markersize', 12, 'LineWidth', 3)

% Find the strike points in the standard EFIT

Rlessthan = limdata(2,:) <  2.0;
Zlessthan = limdata(1,:) < -0.8;
 
boxL = Rlessthan & Zlessthan;
limIdxL = find(boxL ~= 0);

spRZLP = isoflux_spFinder(psizr257257, psixPL, rg257, zg257, limdata, limIdxL);

plot(spRZLP(1,1), spRZLP(1,2), 'om', 'LineWidth', 2)
plot(spRZLP(2,1), spRZLP(2,2), 'om', 'LineWidth', 2)

% Determine major radius of primary separatrix on the midplane
 
psimid = interp2(rg257, zg257, psizr257257, rg257, zmaxis);
 
idxP = find(psimid > psixPL, 1, 'last' ); % midplane primary

idxPrimary   = [idxP-1 idxP idxP+1 idxP+2];

ppPrimary   = spline(rg257(idxPrimary),   psimid(idxPrimary));

cP = ppPrimary.coefs(2,:);

rootsP = roots([cP(1) cP(2) cP(3) cP(4)-psixPL]) + rg257(idxP);

[~,idxminP] = min(abs(rootsP-rg257(idxP)));

rmidPrimary   = rootsP(idxminP);

% Compute the poloidal and total fields at the midplane

BrMid = interp2(rg257, zg257, Br257257, rmidPrimary, zmaxis);
BzMid = interp2(rg257, zg257, Bz257257, rmidPrimary, zmaxis);

BpMid = sqrt(BrMid*BrMid + BzMid*BzMid);

BTMid = (bzero*rzero)/rmidPrimary;

BTotMid = sqrt(BpMid*BpMid + BTMid*BTMid);

% Compute the poloidal and total fields at the outer strike point

BrSP = interp2(rg257, zg257, Br257257, spRZLP(2,1), spRZLP(2,2));
BzSP = interp2(rg257, zg257, Bz257257, spRZLP(2,1), spRZLP(2,2));

BpSP = sqrt(BrSP*BrSP + BzSP*BzSP);

BTSP = (bzero*rzero)/spRZLP(2,1);

BTotSP = sqrt(BpSP*BpSP + BTSP*BTSP);

% Compute the flux expansion at the outer strike point

fExp = (BpMid/BTotMid)/(BpSP/BTotSP);

% Load the IRTV heat flux data

qperp = transpose(importdata('qperp_158139_2850.txt'));

s = importdata('s_158139_2850.txt');

% Index data for outer SP

idxOuter = find(s > 140);
idxInner = setdiff(1:length(s), idxOuter);

% Index data for gap

idxGap = find(s < 170);

% Remove the gap

gap1 = s(idxGap(end));
gap2 = s(idxGap(end)+1);

dgap = gap2 - gap1;

s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;

% Fit Eich profile to the data

qEich = @(q0, S, lambdaQ, s0, qBG, fExp, x) ...
    (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBG;

ft = fittype(qEich, 'problem', 'fExp', 'independent', 'x');

options = fitoptions(ft);
options.StartPoint = [950 1.4 0.2 155 40];
options.Lower = [0 0 0 0 0];
 
fitEich = fit(s(idxOuter), qperp(idxOuter), ft, 'problem', fExp, options);
 
% Extract the coefficients
 
q0      = fitEich.q0;
S       = fitEich.S;
lambdaQ = fitEich.lambdaQ;
s0      = fitEich.s0;
qBG     = fitEich.qBG;

% Define the optimized Eich profile

s_Eich = linspace(0,300,1000);

qperp_Eich = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (s_Eich - s0)./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (s_Eich - s0)./S) + qBG;

% Plot the data and fit

figure(2)
plot(s, qperp, '-ob', 'LineWidth', 1, 'MarkerSize', 2)

hold on
plot(s_Eich, qperp_Eich, '-r', 'LineWidth', 1)

xlim([min(s) max(s)])

xlabel('Distance s [cm]')
ylabel('Heat Flux [W/cm^2]')

title([int2str(158139) ': ' int2str(2850) ' ms'])

% Add fit parameters to plot

text(100, 450, ['S = ' num2str(S,'%3.1f') ' cm'], 'FontSize', 12)

text(100, 400, ['\lambda_{q} = ' num2str(lambdaQ,'%3.1f') ' cm'], ...
    'FontSize', 12)

text(100, 350, ['f_{exp} = ' num2str(fExp,'%3.1f')], 'FontSize', 12)
