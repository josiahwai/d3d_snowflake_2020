function qperp = gaussProfileWithBG(s, q0, s0, sigma, qBG)

qperp = q0*exp(-(s-s0).^2/sigma) + qBG;

end
