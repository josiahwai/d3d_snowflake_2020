close all

eq = eqs{end};
struct_to_ws(eq);

[rxP, rxS, zxP, zxS, psixP, psixS] = my_snowfinder(rg, zg, psizr, psibry);
plot_eq(eq);
axis([1 1.5 -1.6 -.9]);
scatter([rxP rxS], [zxP zxS], 'filled')