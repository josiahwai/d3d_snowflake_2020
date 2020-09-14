ccc

load('psi.mat')
load('q0.mat')
load('rdiv.mat')
load('sdiv.mat')
load('tau.mat')
load('zdiv.mat')
load('thetaB')

D = .2;

heat_diffusion(q0, sdiv, tau, D)

function [qpar, qperp] = heat_diffusion(q0, sdiv, tau, D)

n = length(sdiv);

for i = 1:n
  si = sdiv(i);
  t = tau(i);
  
  % green's function for heat equation
  greens =  1 / sqrt(4*pi*t) * exp( - (sdiv - si).^2 / (4 * D * t));         
  
  q(i) = trapz(sdiv, q0 .* greens);
  
end

% qperp = q.*sin(thetaB);

open('qir.fig')
hold on
plot(sdiv, q * 0.13/ max(q), '--b')    
% plot(sdiv, qperp * 0.13/ max(qperp), 'b')    

end


