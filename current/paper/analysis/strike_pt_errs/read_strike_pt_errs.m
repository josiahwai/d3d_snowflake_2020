load('strike_sfm.mat')


% remove outliers (usually errors in fitting), nans, bad data
[~,i1] = rmoutliers(strike.err0, 'ThresholdFactor', 3 );
i2 = sum(abs(strike.errf), 2) > 0.1;
i3 = sum ( isnan( strike.err0), 2);
i4 = sum ( isnan( strike.errf), 2);

iuse = ~(i1 | i2 | i3 | i4);

% strike point errors
err0 = strike.err0(iuse,:);
errf = strike.errf(iuse,:);

mean(err0)
mean(errf)

sum(mean(err0))
sum(mean(errf))


% power fraction
ipf = logical( sum( isnan( [strike.pf_true' strike.pf_0' strike.pf_f']), 2));
iuse_pf = ~(i1 | i2 | i3 | i4 | ipf);

pf_true = strike.pf_true(iuse_pf);
pf_0 = strike.pf_0(iuse_pf);
pf_f = strike.pf_f(iuse_pf);


mean(abs(pf_0 - pf_true))
mean(abs(pf_f - pf_true))



% RESULTS
%  2% err power fraction vs 16%
%  7.7 cm vs 1.1 cm strike pts (summed error, sfm)
%  5.4cm vs 0.4cm strike pts (summed error, sfp)











