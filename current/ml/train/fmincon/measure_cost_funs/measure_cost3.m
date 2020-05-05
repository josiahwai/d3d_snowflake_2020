function J = measure_cost3(sim)

% ---------------------------------------
  % QUADRATIC COSTS ON MATCH TO HEAT FLUX
  % -------------------------------------
  J = 0.0;  
  struct_to_ws(sim);
  
  wt_pk = [1 1 20];
  wt_q  = [5 5 5] * 0.002;  
  
  % which peaks to include in cost analysis
  usepkI = ~isnan(s_qirmax(1));
  usepkX = (~isnan(s_qirmax(2)) | ~isnan(s_qmax(2))) &...
    sim.qirmaxN(2) > .025;
  usepkO = ~isnan(s_qirmax(3));
  
  ipk = boolean([usepkI usepkX usepkO]);
  
  % convert nans to zero
  qmaxN(isnan(qmaxN)) = 0;
  qirmaxN(isnan(qirmaxN)) = 0;
  s_qmax(isnan(s_qmax)) = 0;
  s_qirmax(isnan(s_qirmax)) = 0;
  
  % cost on peak distances
  J = J + sum(wt_pk(ipk) .* (s_qmax(ipk) - s_qirmax(ipk)).^2);
  
  % cost on peak relative magnitudes
  J = J + sum(wt_q(ipk) .* (qmaxN(ipk) - qirmaxN(ipk)).^2);
  
  
  J = double(J);  
        


