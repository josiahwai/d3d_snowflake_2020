
% rxPL =  1.1202;
% zxPL = -1.0643;
% 
% rxSL =  1.1438;
% zxSL = -1.2796;

rxP_PRED =  1.1302 - 0.01;
zxP_PRED = -1.1142 + 0.05;

rxS_PRED =  1.1838 - 0.04;
zxS_PRED = -1.2996 + 0.02;

for ii = 59:63
    
    zxS_PRED = zxS_PRED + 0.002; 
    
%     rxP_PRED = rxP_PRED + 0.01*cos(pi/4);
%     zxP_PRED = zxP_PRED + 0.01*sin(pi/4);
    
%     rxS_PRED = rxS_PRED - 0.001;
%     zxS_PRED = zxS_PRED + 0.001; 
    
    gsdesign_d3d_165288_4000_EQHF
    
    save(['eq_' int2str(ii) '.mat'], 'eq')
    
    close all
    
end

