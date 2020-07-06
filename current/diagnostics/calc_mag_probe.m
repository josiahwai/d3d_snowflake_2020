
% Calculate B-field experienced by the magnetic probes, for a given
% equilibrium. eq = equilibrium with flux map, xmp=radius of mag probes,
% ymp = height of mag probes, amp = angle orientation of mag probes in deg

% outputs: Bprobe = calculated B-field (ie, what an ideal probe would measure
% if 100% accurate and equlibrium is 100%  accurated), chisq = f

function Bprobe = calc_mag_probe(eq, xmp, ymp,amp)

if isfield(eq, 'gdata'), eq = eq.gdata; end
struct_to_ws(eq);

[~, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, xmp, ymp);

Br = -1./(2*pi*xmp).*psi_z;
Bz =  1./(2*pi*xmp).*psi_r;       
Bprobe = Br.*cosd(amp) + Bz.*sind(amp);









