function [qmax, s_qmax, r_qmax, z_qmax, psi_qmax, psir_qmax, psiz_qmax] = qpeak_info(...
    q, s, pkthresh, sLimTot, limdata, rg, zg, psizr)

try   
    [qmax,k] = findpeaks(q,'NPeaks',1,'sortstr','descend',...
        'minpeakheight', pkthresh, 'minpeakprominence', pkthresh);
    
    if isempty(k)
        [qmax,s_qmax,r_qmax,z_qmax,psi_qmax] = unpack(nan(5,1));
        return
    end
    
    s_qmax = s(k);
    [r_qmax, z_qmax] = calcLimDistanceInv(sLimTot-s_qmax,limdata);
    [psi_qmax, psir_qmax, psiz_qmax] = bicubicHermite(rg, zg, psizr, r_qmax, z_qmax);
catch
end