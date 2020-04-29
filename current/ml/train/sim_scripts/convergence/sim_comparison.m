load('/u/jwai/d3d_snowflake_2020/current/optimize/convergence/active-set/history.mat')
as = history;
load('/u/jwai/d3d_snowflake_2020/current/optimize/convergence/interior-point/history.mat')
ip1 = history;
load('/u/jwai/d3d_snowflake_2020/current/optimize/convergence/interior-point2/history.mat')
ip2 = history;
load('/u/jwai/d3d_snowflake_2020/current/optimize/convergence/history.mat')
ip3 = history;
ip3 = ip2;

set(0,'defaultlinelinewidth',2);
plot(as.fval)
hold on
plot(ip1.fval)
plot(ip2.fval)
plot(ip3.fval)
legend('as','ip1','ip2','ip3')


figure()
load('/u/jwai/d3d_snowflake_2020/current/inputs/tok_data/d3d_obj_mks_struct_6565.mat')
limdata = tok_data_struct.limdata;
plot(limdata(2,:), limdata(1,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
hold on
axis equal
axis([1.0 1.5 -1.4 -0.9])

% plot(as.x(:,1), as.x(:,3),'b')
% plot(as.x(:,2), as.x(:,4),'b')
% 
% plot(ip1.x(:,1), ip1.x(:,3),'r')
% plot(ip1.x(:,2), ip1.x(:,4),'r')
% 
% plot(ip2.x(:,1), ip2.x(:,3),'r')
% plot(ip2.x(:,2), ip2.x(:,4),'r')
% 
% plot(ip3.x(:,1), ip3.x(:,3),'g')
% plot(ip3.x(:,2), ip3.x(:,4),'g')


r = [as.x(end,1:2) ip1.x(end,1:2) ip2.x(end,1:2) ip3.x(end,1:2)];
z = [as.x(end,3:4) ip1.x(end,3:4) ip2.x(end,3:4) ip3.x(end,3:4)];
plot(r, z, 'xk', 'Markersize', 10, 'LineWidth', 4)

% 1 STEP OPTIMIZATION
plot([1.1267    1.1475], [-1.1739   -1.2828], 'xr', 'Markersize', 10, 'LineWidth', 4)

% INITIAL EFIT
plot([1.1056 1.1603], [-1.2072 -1.2597], 'xb', 'Markersize', 10, 'LineWidth', 4)




















