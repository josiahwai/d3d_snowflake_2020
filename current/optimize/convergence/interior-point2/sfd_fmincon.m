function [history,searchdir] = sfd_fmincon

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
clear designeq_cost

% call optimization 
x0 = [1.1267    1.1472   -1.1665   -1.2909]; % predicted optimum from dflux1

lb = x0-0.1;
ub = x0+0.1;

options = optimoptions(@fmincon,'OutputFcn',@outfun,...
    'Display','iter','Algorithm','interior-point');

xsol = fmincon(@designeq_cost,x0,[],[],[],[],lb,ub,[],options);
save('history','history')


x = history.x;
plot(x(:,1),x(:,3),'k')
plot(x(:,2),x(:,4),'b')


function stop = outfun(x,optimValues,state)
    stop = false;

    switch state
        case 'init'            
        case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x];
        case 'done'
        otherwise
    end
end
end











