function J = measure_cost2(sim)

J = 0.0;
struct_to_ws(sim);


% resample sir, qir for equally spaced sir
% .......................................
[sir_old,k] = unique(sir);
sir = linspace(sir_old(1), sir_old(end), length(sir_old));
qir = interp1(sir_old,qir(k),sir);
qir = qir/max(qir);

% resample s,q to be same size as sir,qir
% .......................................

% remove duplicates
[sI,k] = unique(sI); qI = qI(k);
[sX,k] = unique(sX); qX = qX(k);
[sO,k] = unique(sO); qO = qO(k);

s = sir;
q = zeros(1,length(s));

% resample
for k = 1:length(s)
  if s(k) > min(sI) && s(k) < max(sI)
    q(k) = interp1(sI,qI,s(k));
  elseif s(k) > min(sX) && s(k) < max(sX)
    q(k) = interp1(sX,qX,s(k));
  elseif s(k) > min(sO) && s(k) < max(sO)
    q(k) = interp1(sO,qO,s(k));
  end
end

q = q/max(q);
q(q==0) = median(qir);

% close all
% plot(s,q)
% hold on
% plot(sir,qir)


[~,i_qmax]   = min(abs(s'-s_qmax));
[~,i_qirmax] = min(abs(sir'-s_qirmax));

dqA = q(i_qmax) - qir(i_qmax);
dqB = q(i_qirmax) - qir(i_qirmax);

wt3 = [1 1 10];
wt2 = wt3([1 3]);

J = 0;
if length(dqA) == 3
  J = J + sum(wt3.*dqA.^2);
else
  J = J + sum(wt2.*dqA.^2);
end

if length(dqA) == 3
  J = J + sum(wt3.*dqB.^2);
else
  J = J + sum(wt2.*dqB.^2);
end

J = double(J);











