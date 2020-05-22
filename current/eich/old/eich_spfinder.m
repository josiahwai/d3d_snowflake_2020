% finds the strike points from the heat flux data

function [rsp,zsp,ssp] = eich_spfinder(s,q,tok_data_struct,plotit)

eichfun_i = @(q0, S, lambdaQ, fExp, s0, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (-(x - s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x - s0))./S);

eichfun_o = @(q0, S, lambdaQ, fExp, s0, qBg, x) (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
  - (x-s0)./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (x-s0)./S) + qBg;

if size(s,2) ~= 1, s = s'; end
ii = find(s<120);
io = find(s>145);
ix = find(s>=120 & s<=145);

% fit the outer peak to find strike pt
% ....................................

% remove gap before fiting
iGap = find(s < 170,1,'last');
dgap = s(iGap+1) - s(iGap);
s_nogap = s;
s_nogap(iGap(end)+1:end) = s_nogap(iGap(end)+1:end) - dgap;


ft = fittype(eichfun_o);
options = fitoptions(ft);
options.StartPoint = [100 1 1 4 150 4];
options.Lower = [0 0 0 0 0 0];  

fito = fit(s_nogap(io), q(io), ft, options);
sspo =  fito.s0;


% fit the middle peak, if the peak is there
% .........................................
pkthresh = 1.5*median(q);

[qpkx,ipkx] = findpeaks(q(ix), 'NPeaks',1,'sortstr','descend',...
  'minpeakheight', pkthresh, 'minpeakprominence', pkthresh);

% estimate strike point as location where q = 3/4*q_peak on inbd side
[~,k] = min(abs(q(ix(1:ipkx)) - 3/4*qpkx));
sspx = s(ix(k));

if isempty(sspx), sspx = nan; end

% fit the inner peak to find strike pt
% ....................................
ft = fittype(eichfun_i);
options = fitoptions(ft);
options.StartPoint = [100 1 1 4 150];
options.Lower = [0 0 0 0 0];  

fiti = fit(s(ii), q(ii), ft, options);
sspi =  fiti.s0;

% convert from s to (r,z)
ssp = [sspi sspx sspo]/100;

limdata = tok_data_struct.limdata;
slimtot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);
[rsp,zsp] = calcLimDistanceInv(slimtot - ssp, limdata);

if plotit
  figure
  hold on
  plot(s,q)
  xline(sspi);
  if ~isnan(sspx), xline(sspx); end
  xline(sspo);
end
end









