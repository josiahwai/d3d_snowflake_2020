function q = eichfun_inner(S,fExp,lambdaQ,q0,s0,x)

q = (q0/2)*exp((S/(2*lambdaQ*fExp))^2 ...
      - (-(x - s0))./(lambdaQ*fExp)).*erfc(S/(2*lambdaQ*fExp) - (-(x - s0))./S);
end

    
