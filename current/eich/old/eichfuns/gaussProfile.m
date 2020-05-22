function qperp = gaussProfile(s, q0, s0, sigma)

qperp = q0*exp(-(s-s0).^2/sigma);

end
