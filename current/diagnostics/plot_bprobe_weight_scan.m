close all

wts = logspace(-2,6,30);
wts = [nan wts(3:end)];

load('bprobe_weight_scan_data.mat')

% xp targets
rx0 = [1.2235    1.1560];
zx0 = [-1.1746   -1.3229];
load('xps');
xpf = xps{end};
rxf =xpf(1:2);
zxf = xpf(3:4);

% psi targets
load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm/155354_sfm/3727/eqs.mat')
snow_tf = analyzeSnowflake(eqs{end});
snow_t0 = analyzeSnowflake(eqs{1});

psixPL_tf = snow_tf.psixPL;
psixSL_tf = snow_tf.psixSL;
dpsi_tf = psixPL_tf - psixSL_tf;

psixPL_t0 = snow_t0.psixPL;
psixSL_t0 = snow_t0.psixSL;
dpsi_t0 = psixPL_t0 - psixSL_t0;


% from eich_fit to 155354:3727
ef.rsp = [1.0160 1.0908 1.2572];
ef.zsp = [-1.0756 -1.2995 -1.3630];


% LOAD DATA
      
dxp = [];
dxp_targ = [];
is_sfp = [];
dsp = [];
dpsi_err = [];
psixSL = [];
psixPL = [];

plot_eq(bprobe_weight_scan_data(1).eqs)

for k = 1:length(bprobe_weight_scan_data)
  ssq(k) = bprobe_weight_scan_data(k).ssq;
  eq = bprobe_weight_scan_data(k).eqs;
  snow = analyzeSnowflake(eq);
  
  % x-points
  rx = snow.rx;
  zx = snow.zx;  
  dxp(k) = sum (sqrt ( (rx-rx0).^2 + (zx-zx0).^2)); 
  dxp_targ(k) = sum (sqrt ( (rx-rxf).^2 + (zx-zxf).^2));     
  % scatter(snow.rx,snow.zx,'filled')    
  
  
  % delta flux
  dpsi = snow.psixPL - snow.psixSL;
  dpsi_err(k) = abs(dpsi_tf - dpsi);  
  psixPL(k) = snow.psixPL;
  psixSL(k) = snow.psixSL;
  
         
  % strike point targets
  isp = [snow.rSPP(1) snow.zSPP(1)];
  
  if snow.snowPlus
    osp = [snow.rSPP(end) snow.zSPP(end)]; 
    xsp = [snow.rSPS(2) snow.zSPS(2)];
  else
    osp = [snow.rSPS(end) snow.zSPS(end)];
    xsp = [snow.rSPP(2) snow.zSPP(2)];
  end
      
  rsp = [isp(1) xsp(1) osp(1)];
  zsp = [isp(2) xsp(2) osp(2)];
  
  dsp(k) = sum (sqrt ( (rsp-ef.rsp).^2 + (zsp-ef.zsp).^2)); 
        
  % snow type
  is_sfp(k) = snow.snowPlus;
  
  
end

%%
close 

figure
hold on
wt = ((1:length(ssq)) - 15)/10;
wt(wt<0) = nan;

wt_sfp = (20.5 - 15) / 10; % transition from sfp to sfm 

colorblind_cmap


% yyaxis left
plot(wt, dxp_targ*100 , 'color', blugreen, 'linewidth', 1.5);  % [cm]
plot(wt, dsp*100, 'color', orange, 'linewidth', 1.5 );  % [cm]
plot(wt, dpsi_err*1000, 'color', skyblue, 'linewidth', 1.5 ); % [mWB]
plot(wt, sqrt(ssq)*100, 'color', redpurple, 'linewidth', 1.5)  % [100*Gs]

plot([.2615 .3], [9.2137 9.85], 'color', black, 'linewidth', 0.75)

xline(wt_sfp, 'color', gray, 'linewidth', 1);

xlim([0 1])
ylim([ 0 12.5])
box('on')

t1 = text(0.01, 0.68, '\Deltaxp [cm]', 'units', 'normalized', 'color', blugreen);
t2 = text(0.01, 0.91, '\Deltasp [cm]', 'units', 'normalized', 'color', orange);
t3 = text(0.28, 0.82, '\Delta\psi [mWb]', 'units', 'normalized', 'color', skyblue);
t4 = text(0.01, 0.23, '\DeltaBprobe [100Gs]', 'units', 'normalized', 'color', redpurple);

set([t1 t2 t3 t4], 'fontweight', 'bold', 'fontsize', 12)


ylabel('\Sigma errors', 'fontsize', 12, 'fontweight', 'bold')
xlabel('IRTV weight', 'fontsize', 12, 'fontweight', 'bold')

set(gcf, 'position', [881 280 560 420])

x1 = [wt_sfp wt_sfp - 0.2];
y = [1 1] * 11.3;

drawArrow(x1,y,{'Color',gray,'LineWidth',4}); 
t = text(wt_sfp - 0.025, 11.85, 'SFD+', 'fontweight', 'bold', 'fontsize', 14, ...
  'color', gray, 'horizontalAlignment', 'right');

t = text(wt_sfp + 0.025, 11.85, 'SFD-', 'fontweight', 'bold', 'fontsize', 14, ...
  'color', gray, 'horizontalAlignment', 'left');
x2 = [wt_sfp wt_sfp + 0.2];
drawArrow(x2,y,{'Color',gray,'LineWidth',4}); 

title('155354: 3720ms', 'fontsize', 14)
fn = '/u/jwai/d3d_snowflake_2020/current/diagnostics/fig_bprobes.eps';
saveas(gcf, fn, 'epsc')




