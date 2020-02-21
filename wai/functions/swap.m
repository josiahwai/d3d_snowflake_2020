function swap(A,B)
% SWAP(A,B) puts the contents of variable A into variable B and vice versa.

narginchk(2,2)

evalstr = sprintf('[%s,%s] = deal(%s,%s) ;',inputname(2),inputname(1),inputname(1),inputname(2));
evalin('caller', evalstr)



