function add_lim_colorcode(fignum)


rv = [1.016 1.016 1.153 1.42 1.42 1.372 1.372 1.768];
zv = [-0.83 -1.223 -1.363 -1.363 -1.329 -1.329 -1.25 -1.25 -1.25];
nv = length(rv);
cmap = colormap(lines);

% figure(19)
% hold on
% for i = 1:nv-1
%   plot( rv(i:i+1), zv(i:i+1), 'color', cmap(i,:), 'linewidth', 4)
% end
% set(gcf,'position', [1569 607 350 271])

   
load('d3d_obj_mks_struct_6565.mat')
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

sv = sLimTot - calcLimDistance(rv, zv, limdata);
sv = sv*100;

figure(fignum)
hold on
for i = 1:nv-1
   plot( sv(i:i+1), [0 0], 'color', cmap(i,:), 'linewidth', 3)
end
   
   
   
   
   