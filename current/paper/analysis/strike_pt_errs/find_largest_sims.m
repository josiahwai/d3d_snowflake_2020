e = sum(strike.err0 - strike.errf,2);
[dist, i] = maxk(e, 10)

[strike.shot(i)' strike.time(i)']



155348        3760


