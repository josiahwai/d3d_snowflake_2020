clear; clc; close all

wts = logspace(1,4.7,100);


% UNCOMMENT FOR 155354:3727
load('bprobe_weight_scan_data.mat')
load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm_large_lambdaq/155354_sfm/3727/xps.mat')
load('/u/jwai/d3d_snowflake_2020/current/sfmodel/jobs/sfm_large_lambdaq/155354_sfm/3727/eqs.mat')
% from eich_fit to 155354:3727
ef.rsp = [1.0160 1.0908 1.2742];
ef.zsp = [-1.0758 -1.2995 -1.3630];

% UNCOMMENT FOR 155350:3893
% load('bprobe_weight_scan_data155350.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155350_3893/xps.mat')
% load('/u/jwai/d3d_snowflake_2020/current/paper/fig_scripts/fig_path/155350_3893/eqs.mat')
% ef.rsp = [1.0160    1.0862    1.2363];
% ef.zsp = [-1.0637   -1.2947   -1.3630];


neqs = length(bprobe_weight_scan_data);
wts  = wts(1:neqs);

% xp targets
rx0 = [1.2235    1.1560];
zx0 = [-1.1746   -1.3229];
xpf = xps{end};
rxf =xpf(1:2);
zxf = xpf(3:4);

% psi targets
snow_tf = analyzeSnowflake(eqs{end});
snow_t0 = analyzeSnowflake(eqs{1});

psixPL_tf = snow_tf.psixPL;
psixSL_tf = snow_tf.psixSL;
dpsi_tf = psixPL_tf - psixSL_tf;

psixPL_t0 = snow_t0.psixPL;
psixSL_t0 = snow_t0.psixSL;
dpsi_t0 = psixPL_t0 - psixSL_t0;

% LOAD DATA
      
dxp = [];
dxp_targ = [];
is_sfp = [];
dsp = [];
dpsi_err = [];
psixSL = [];
psixPL = [];

plot_eq(bprobe_weight_scan_data(1).eqs)

for k = 1:neqs
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
close all
figure
hold on

% normalize the x-axiswt
x0 = 40;
xf = 89;
wt = ((1:neqs) - x0) / (xf - x0);
wt( wt < 0) = nan;
wt( wt > 1) = nan;

is_sfp(1) = 1;
i = find(is_sfp == 0, 1);
wt_sfp = (wt(i) + wt(i-1)) / 2;

colorblind_cmap


% yyaxis left
plot(wt, dxp_targ*100 , 'color', blugreen, 'linewidth', 1.5);  % [cm]
plot(wt, dsp*100, 'color', orange, 'linewidth', 1.5 );  % [cm]
plot(wt, dpsi_err*1000, 'color', skyblue, 'linewidth', 1.5 ); % [mWB]
plot(wt, sqrt(ssq)*200, 'color', redpurple, 'linewidth', 1.5)  % [50*Gs]

% plot([.2615 .3], [9.2137 9.85], 'color', black, 'linewidth', 0.75)

xline(wt_sfp, 'color', gray, 'linewidth', 1);

% xlim([0 1])
ylim([ 0 15])
box('on')
grid on

t1 = text(0.01, 0.76, '\Deltaxp [cm]', 'units', 'normalized', 'color', blugreen);
t2 = text(0.01, 0.52, '\Deltasp [cm]', 'units', 'normalized', 'color', orange);
t3 = text(0.01, 0.86, '\Delta\psi [mWb]', 'units', 'normalized', 'color', skyblue);
t4 = text(0.01, 0.25, '\DeltaBprobe [50Gs]', 'units', 'normalized', 'color', redpurple);

set([t1 t2 t3 t4], 'fontweight', 'bold', 'fontsize', 12)


ylabel('\Sigma errors', 'fontsize', 12, 'fontweight', 'bold')
xlabel('IRTV weight', 'fontsize', 12, 'fontweight', 'bold')

set(gcf, 'position', [211 254 560 420])

x1 = [wt_sfp wt_sfp - 0.2];
y = [1 1] * 14.2;

drawArrow(x1,y,{'Color',gray,'LineWidth',4}); 
t = text(wt_sfp - 0.025, 13.3, 'SFD+', 'fontweight', 'bold', 'fontsize', 14, ...
  'color', gray, 'horizontalAlignment', 'right');

t = text(wt_sfp + 0.025, 13.3, 'SFD-', 'fontweight', 'bold', 'fontsize', 14, ...
  'color', gray, 'horizontalAlignment', 'left');
x2 = [wt_sfp wt_sfp + 0.2];
drawArrow(x2,y,{'Color',gray,'LineWidth',4}); 

title('155354: 3720ms', 'fontsize', 14)
fn = '/u/jwai/d3d_snowflake_2020/current/diagnostics/fig_bprobes.eps';
saveas(gcf, fn, 'epsc')




