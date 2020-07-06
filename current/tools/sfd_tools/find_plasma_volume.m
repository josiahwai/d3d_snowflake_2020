% For arbitrary cross-section, the volume of a revolved figure with
% cross-sectional area A, and geometric centroid Rc, is V = 2*pi*Rc*A
% source: https://en.wikipedia.org/wiki/Pappus's_centroid_theorem

function [V, A, Rc] = find_plasma_volume(eq, plotit)

if nargin < 2, plotit = 0; end

struct_to_ws(eq);

bry = polyshape(rbbbs, zbbbs);

[Rc, Zc] = centroid(bry);   % [m]
A = polyarea(rbbbs, zbbbs); % [m^2]
V = 2*pi*Rc*A;              % [m^3]


if plotit
  figure
  hold on
  plot(bry)
  scatter(Rc, Zc, 'filled')
end



