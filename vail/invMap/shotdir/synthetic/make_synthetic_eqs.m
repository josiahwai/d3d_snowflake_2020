
shot = 165288;
time = 4000;

topdir = '/u/pvail/d3d_snowflake_2019/heatfluxconstrain/165288/4000/';

efit_dirname = [topdir 'EFIT01'];
    
init = read_eq(shot, time/1000, efit_dirname);

psizr  = init.gdata.psizr;
psibry = init.gdata.psibry;

rg = init.gdata.rg;
zg = init.gdata.zg;

rExp   =  1.1500;
zExp   = -1.2500;
rhoExp =  0.1000;

[rxPL, zxPL, rxSL, zxSL, ~, ~, ~, ~, ~, ~] = ...
    snowFinder(psizr, rExp, zExp, rhoExp, rg, zg);

[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    
    temp   = psixPL;
    psixPL = psixSL;
    psixSL = temp;
    
    temp1 = rxPL;
    temp2 = zxPL;
    rxPL  = rxSL;
    zxPL  = zxSL;
    rxSL  = temp1;
    zxSL  = temp2;
    
end

for ii = 31:35
    
    rxPL = rxPL + 0.01*cos(pi/4);
    zxPL = zxPL + 0.01*sin(pi/4);
    
    rxSL = rxSL - 0.01*cos(pi/4);
    zxSL = zxSL - 0.01*sin(pi/4);
    
    gsdesign_d3d_165288_4000_EQHF
    
    save(['eq_' int2str(ii) '.mat'], 'eq')
    
    close all
    
end

