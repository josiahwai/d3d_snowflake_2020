

for k = 1:length(eqs)
  eq1 = eqs{k};
  snow1 = analyzeSnowflake(eq1);
  
  ssp = [snow1.sSPP(1:2) snow1.sSPS(end)];
  ef.ssp = [1.0755    1.3300    1.5228];
  
  k
  ssp - ef.ssp
end

1.0800    1.3335    1.5179








