clear; clc; close all;

shot = 155338;

qperp = load(['qperp_' num2str(shot) '.mat']);
s = qperp.s;
t = qperp.t;
qperp = qperp.qperp;

gif_fn = ['/p/nstxusr/nstx-users/jwai/d3d_snowflake_2020/current/inputs/extras/' num2str(shot) '.gif'];

eqdir = '/p/nstxusr/nstx-users/jwai/d3d_snowflake_2020/current/inputs/eqs/efit01/';

d = dir([eqdir num2str(shot)]);
d(1:2) = [];

levels = linspace(.02,-0.02,20);

eq = read_gfile_tok(d(1).name, 'd3d');  

figure
set(gcf, 'Position', [406 501 574 723]);

ax(1) = subplot(3,1,1:2);
hold on
plot(eq.xlim, eq.ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
axis equal
axis([1.0 1.5 -1.45 -0.95])
xlabel('R [m]')
ylabel('Z [m]')

ax(2) = subplot(3,1,3);
grid on
hold on
ylabel('q_\perp [W/cm^2]')
xlabel('s [cm]')

first_pass = true;

for i = 1:length(d)    
    
  eq_t = str2double(d(i).name(end-3:end));
  
  if eq_t > 3500 && eq_t < 3700
        
    eq = read_gfile_tok([d(i).folder '/' d(i).name], 'd3d');           
    
    axes(ax(1))
    if exist('Hc1', 'var')
      delete(Hc1)
      delete(Hc2)
    end
    [~, Hc1] = contour(eq.rg, eq.zg, eq.psizr, eq.psibry + levels, 'k'); 
    [~, Hc2] = contour(eq.rg, eq.zg, eq.psizr, [eq.psibry eq.psibry], 'r', 'linewidth', 2);
    title(d(i).name)


    axes(ax(2)) 
    if exist('H3', 'var')
      delete(H3)    
    end  
    [~,j] = min(abs(eq_t - t));
    H3 = plot(s, qperp(j,:), 'k');
    ylim([0 max(qperp(:))])
    title(d(i).name)
    drawnow
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    
    if first_pass
      first_pass = false;
      imwrite(imind,cm,gif_fn,'gif', 'Loopcount',inf); 
    else
      imwrite(imind,cm,gif_fn,'gif','WriteMode','append'); 
    end
    
    % pause(0.3)
  end
end













