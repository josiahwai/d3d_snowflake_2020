


rv = [1.016 1.016 1.153 1.42 1.42 1.372 1.372 1.768];
zv = [-0.83 -1.223 -1.363 -1.363 -1.329 -1.329 -1.25 -1.25 -1.25];
nv = length(rv);
cmap = colormap(lines);

figure
hold on
for i = 1:nv-1
  plot( rv(i:i+1), zv(i:i+1), 'color', cmap(i,:), 'linewidth', 2)
end

   
load('d3d_obj_mks_struct_6565.mat')
limdata = tok_data_struct.limdata;
sLimTot = calcLimDistance(limdata(2,1), limdata(1,1), limdata);

sv = sLimTot - calcLimDistance(rv, zv, limdata);

figure
hold on
for i = 1:nv-1
   plot( sv(i:i+1), [0 0], 'color', cmap(i,:), 'linewidth', 2)
end
   
   
   
   
   