function [qmax, s_qmax, r_qmax, z_qmax, psi_qmax] = qpeak_info(...
    q,pkthresh,sdiv,rdiv,zdiv,rg,zg,psizr)

% initialize
[qmax,s_qmax,r_qmax,z_qmax,psi_qmax] = unpack(nan(5,1));

[qmax,k] = findpeaks(q,'NPeaks',1,'sortstr','descend',...
            'minpeakheight', pkthresh, 'minpeakprominence', pkthresh);  

[s_qmax,r_qmax, z_qmax] = unpack([sdiv(k) rdiv(k) zdiv(k)]);

psi_qmax = bicubicHermite(rg, zg, psizr, r_qmax, z_qmax);
