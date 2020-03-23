ccc
shot = 165288;
time_ms = 4200;
wt_pk = [1 1 20];
wt_q  = [1 1 10] * 0.002;
wt_xp = 0.001;


root = '/u/jwai/d3d_snowflake_2020/current/';
simdir = [root 'optimize/output/'];
load('xp_list')

% load and analyze original efit_eq
efit_dir = [root 'inputs/eqs/efit01/' num2str(shot)];
eq = read_eq(shot, time_ms/1000, efit_dir);

% load ir heat flux
qperp_dir  = [root 'inputs/qperp/' num2str(shot) '/'];
qperp_data = ['qperp_' num2str(shot) '.mat'];

load([qperp_dir qperp_data])  % loads q, s, and t

[~,k] = min(abs(t-time_ms));
qir = qperp(k,:)'/100;
noisefloor = median(qir);


% find the snowflake
psizr  = eq.gdata.psizr;
psibry = eq.gdata.psibry;
[psizr, rg, zg] = regrid(eq.gdata.rg, eq.gdata.zg, psizr, 257, 257);
[rxPL, zxPL, rxSL, zxSL] = snowFinder(psizr, 1.15, -1.25, 0.1, rg, zg);
[rxPL, zxPL, psixPL] = isoflux_xpFinder(psizr, rxPL, zxPL, rg, zg);
[rxSL, zxSL, psixSL] = isoflux_xpFinder(psizr, rxSL, zxSL, rg, zg);

if abs(psixSL - psibry) < abs(psixPL - psibry)
    swap(psixPL, psixSL);
    swap(rxPL, rxSL);
    swap(zxPL, zxSL);   
end


% Evaluate the cost function for each simulation
J = inf(length(xp_list),1);
sims = {};

for k = 1:length(xp_list)
try    
    fn = ['sim' num2str(k)];
    load([simdir fn]);
    sims{k} = sim;            
    
    % snowflake type based on ir
    snowplus = 0;    
    if isnan(sim.qirmax(2)), snowplus = 1; end
        
    % snowflake type from sim / flux
    snowplus_pred = 0;
    if sim.psix(1) < sim.psix(2), snowplus_pred = 1; end
        
    
    % only evaluate cost if eq is approximately the right type of snowflake
    matchesSnowType = 0;    
    if snowplus == snowplus_pred
        matchesSnowType = 1; 
    elseif snowplus
        % sim predicts snowminus, but the x-peak is small
        if isnan(sim.qmax(2)) || sim.qmax(2) < 1.5*noisefloor 
            matchesSnowType = 1; 
        end
    elseif ~snowplus
        % sim predicts snowplus, but the IR x-peak is small
        if isnan(sim.qirmax(2)) || sim.qirmax(2) < 1.5*noisefloor
            matchesSnowType = 1; 
        end        
    end    
    
    
    if matchesSnowType        
        % inner peak distance
        j = 0;         
        
        % what to include in cost analysis
        usepkI = ~isnan(sim.sir(1)); 
        usepkX = (~isnan(sim.sir(2)) | ~isnan(sim.s(2))) &...
            sim.qirmax(2) > 2*noisefloor;
        usepkO = ~isnan(sim.sir(3));
                
        
        % rescale q peak magntiudes       
        qmax = sim.qmax;
        qmax(isnan(qmax)) = 0;
        qmax = qmax / sum(qmax);
        
        qirmax = sim.qirmax;
        qirmax(isnan(qirmax)) = 0;
        qirmax = qirmax / sum(qirmax);
        
        ipk = boolean([usepkI usepkX usepkO]);
        
        % peak distances
        j = j + sum(wt_pk(ipk) .* (sim.s(ipk) - sim.sir(ipk)).^2);
        
        % peak relative magnitudes
        j = j + sum(wt_q(ipk) .* (qmax(ipk) - qirmax(ipk)).^2);
          
        % x-pt movement
        dr = [sim.rx(1)-rxPL sim.rx(2)-rxSL];
        dz = [sim.zx(1)-zxPL sim.zx(2)-zxSL];    
        j = j + wt_xp * sum(dr.^2 + dz.^2);
        
        J(k) = j;
    end
catch
end
end


[cost,iBest] = sort(J);
iBest = iBest(cost~=inf);

k = iBest(1:10);
rxP = xp_list(k,1); 
rxS = xp_list(k,2);
zxP = xp_list(k,3); 
zxS = xp_list(k,4);

rstar = [mean(rxP) mean(rxS)];
zstar = [mean(zxP) mean(zxS)];

openfig('zones.fig');
scatter(rxP,zxP,'r','filled')
scatter(rxS,zxS,'r','filled')
scatter(rstar,zstar,200,'bp','filled')






















