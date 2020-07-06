function regression(param, djmax, rmout, center)

close all
if size(djmax,2) ~= 1, djmax = djmax'; end
if size(param,1) ~= length(djmax), param = param'; end

if rmout
  [~,iout1] = rmoutliers(djmax);
  [~,iout2] = rmoutliers(param);
  iuse = ~(iout1 | iout2);
else
  iuse = true(size(djmax));
end

if center
  param = param(iuse) - mean(param(iuse));
  djmax = djmax(iuse) - mean(djmax(iuse));
else
  param = param(iuse);
  djmax = djmax(iuse);
end

w = pinv(param)*djmax;

dj_pred = param*w;

figure
subplot(1,2,1)
hold on
scatter(djmax, dj_pred)
% axis equal
axis([-2 2 -2 2]*1e5)
xlabel('dj true')
ylabel('dj pred')


subplot(1,2,2)
hold on
dz = linspace(min(param), max(param), 100);
dj = w*dz;
scatter(param, djmax)
plot(dz,dj, '--k')
xlabel('dj true')
ylabel('param')

set(gcf,'position', [440 442 592 261])









