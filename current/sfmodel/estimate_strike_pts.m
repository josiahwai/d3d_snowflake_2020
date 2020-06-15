% only implemented for snowflake plus 

function sp1 = estimate_strike_pts(eq0,sim0)

load('d3d_obj_mks_struct_6565.mat')
limdata = tok_data_struct.limdata;
slim = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

snow0 = analyzeSnowflake(eq0);

if snow0.snowPlus  
  ssp = snow0.sSPP + sim0.s_qirmax([1 3]) - sim0.s_qmax([1 3]);  
  for k = 1:2
    [rsp(k), zsp(k)] = calcLimDistanceInv(slim - ssp(k), limdata);
  end
end

sp1 = [rsp zsp];  















