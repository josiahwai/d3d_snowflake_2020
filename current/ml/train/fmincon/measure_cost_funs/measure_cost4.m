function J = measure_cost4(sim)

% ---------------------------------------
  % QUADRATIC COSTS ON MATCH TO HEAT FLUX
  % -------------------------------------
  J = 0.0;  
  struct_to_ws(sim);
  
  wt_s = [1 1 20];
  wt_q  = [1 1 10] * 0.002;  
  
  % inner peak cost
  % ...............
  
  % ir data has valid inner peak
  if ~isnan(s_qirmax(1)) 
    if ~isnan(s_qmax(1))   % sim peak valid
      J = J + wt_s(1) * (s_qmax(1) - s_qirmax(1)).^2;
      J = J + wt_q(1) * (qmaxN(1) - qirmaxN(1)).^2;
    else 
      J = J + wt_s(1)*.10^2;       % invalid sim peak, 10cm penalty
      J = J + wt_q(1) * qirmaxN(1).^2;
    end
  end
  
  
  
  % outer peak cost
  % ...............
  if ~isnan(s_qirmax(3)) 
    if ~isnan(s_qmax(3))         % sim peak valid
      J = J + wt_s(3) * (s_qmax(3) - s_qirmax(3)).^2;
      J = J + wt_q(3) * (qmaxN(3) - qirmaxN(3)).^2;
    else 
      J = J + wt_s(3)*.10^2;       % invalid sim peak, 10cm penalty
      J = J + wt_q(3) * qirmaxN(3).^2;
    end
  end
  
  
  
  % Middle (X) peak cost
  % ...................
  valid_xpk_ir = ~isnan(s_qirmax(2));
  valid_xpk_sim = ~isnan(s_qmax(2));
  
  % sim and ir have middle peaks. Measure cost as normal
  if valid_xpk_ir && valid_xpk_sim  
    J = J + wt_s(2) * (s_qmax(2) - s_qirmax(2)).^2;
    J = J + wt_q(2) * (qmaxN(2) - qirmaxN(2)).^2;
  % mismatch: measure cost with a 10 cm penalty on position
  elseif valid_xpk_ir ~= valid_xpk_sim 
    J = J + wt_s(2) * .10^2;
    J = J + wt_q(2) * (nan2zero(qmaxN(2)) - nan2zero(qirmaxN(2))).^2;
  end
  
  
  
  if isnan(J)
    J = 100;
    warning('Cost function evaluation resulted in NaN')
  end
  
  J = double(J);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
        


