% Analyze if the new equilibrium was better in predicting heatflux peak
% locations

ccc
plotit = 1;
times = 1000:100:4900;
shot = 155355;


root = '/u/jwai/d3d_snowflake_2020/current/';
addpath(genpath(root));
warning('off','all')

% load ir heat flux for shot
irdir = [root 'inputs/qperp/' num2str(shot) '/'];
fn_ir = ['qperp_' num2str(shot) '.mat'];
load([irdir fn_ir]) % loads qperp, t

sIR = [];
sSim_c = [];
sSim_u = [];
sSim_e = [];

for iTime = 1:length(times)
    try       
        time_ms = times(iTime);
        
        [~,k] = min(abs(t-time_ms));
        q = qperp(k,:);
        
        % Remove the gap from s (distance along limiter)
        idxGap = find(s < 170);
        gap1 = s(idxGap(end));
        gap2 = s(idxGap(end)+1);
        dgap = gap2 - gap1;
        s(idxGap(end)+1:end) = s(idxGap(end)+1:end) - dgap;
        
        iI = find(s<115);
        iO = find(s>145);
        iX = setdiff(1:length(s), [iI iO]);
        
        % find inner peak
        [~,k] = findpeaks(q(iI),'NPeaks',1,'sortstr','descend','minpeakdistance',...
            10,'minpeakheight',0.05*max(q), 'minpeakprominence', 0.05*max(q));
        if ~isempty(k)
            ipkI = min(iI)-1+k;
            sI = s(ipkI);
        else
            sI = 0;
        end
        
        % find region between x-pts peak
        [~,k] = findpeaks(q(iX),'NPeaks',1,'sortstr','descend','minpeakdistance',...
            10,'minpeakheight',0.05*max(q), 'minpeakprominence', 0.05*max(q));
        ipkX = min(iX)-1+k;
        if ~isempty(k)
            ipkX = min(iX)-1+k;
            sX = s(ipkX);
        else
            sX = 0;
        end
        
        % find outboard peak
        [~,k] = findpeaks(q(iO),'NPeaks',1,'sortstr','descend','minpeakdistance',...
            10,'minpeakheight',0.05*max(q), 'minpeakprominence', 0.05*max(q));
        ipkO = min(iO)-1+k;
        if ~isempty(k)
            ipkO = min(iO)-1+k;
            sO = s(ipkO);
        else
            sO = 0;
        end
        
        % plot
        if time_ms == 2600
            figure(1)
            plot(s,q)
            hold on
            ipks = [ipkI ipkX ipkO];
            scatter(s(ipks),q(ipks),'r','filled')
        end
        
        % load heat flux sim data for cake(c), cake unconstrained(u), efit(e)
        hfdir = [root 'outputs/hfsims/'];
        fn_c = ['cake_constrained/' num2str(shot) '/hfsim_' num2str(shot) '_' num2str(time_ms) '.mat'];
        fn_u = ['cake_unconstrained/' num2str(shot) '/hfsim_' num2str(shot) '_' num2str(time_ms) '.mat'];
        fn_e = ['efit_unconstrained/' num2str(shot) '/hfsim_' num2str(shot) '_' num2str(time_ms) '.mat'];
        
        c = load([hfdir fn_c]);
        u = load([hfdir fn_u]);
        e = load([hfdir fn_e]);
        
        % save hf peak locations
        sIR    = [sIR sI sX sO];
        sSim_c = [sSim_c 100*[c.hfsim.sSPI c.hfsim.sSPX c.hfsim.sSPO]];
        sSim_u = [sSim_u 100*[u.hfsim.sSPI u.hfsim.sSPX u.hfsim.sSPO]];
        sSim_e = [sSim_e 100*[e.hfsim.sSPI e.hfsim.sSPX e.hfsim.sSPO]];
    catch
    end
end

% convert nans to zero
sSim_c(isnan(sSim_c)) = 0;
sSim_u(isnan(sSim_u)) = 0;
sSim_e(isnan(sSim_e)) = 0;

figure(2)
hold on
scatter(sIR, sSim_c, 45, 'r', 'o', 'filled')
% scatter(sIR, sSim_u, 45, 'g', 'd', 'filled')
scatter(sIR, sSim_e, 45, 'b', '<', 'filled')
plot(s,s,'--k')

axis([80 200 80 200])
legend('cake constrained', 'efit', 'location', 'northwest')





















