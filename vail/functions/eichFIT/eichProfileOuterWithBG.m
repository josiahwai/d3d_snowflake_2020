function qperp = eichProfileOuterWithBG(s, q0, S, lambdaQ, fExp, s0, qBG)

qperp = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 - (s - s0)./...
    (lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (s - s0)./S) + qBG;

end
