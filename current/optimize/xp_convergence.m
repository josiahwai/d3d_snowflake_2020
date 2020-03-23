% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
ccc
fun = @(x)1+x(1)./(1+x(2)) - 3*x(1).*x(2) + x(2).*(1+x(1));
lb = [0,0];
ub = [1,2];
A = [];
b = [];
Aeq = [];
beq = [];
x0 = (lb + ub)/2;
options = optimoptions(@fmincon,'Display','iter-detailed');

opts=optimoptions(@fmincon,'OutputFcn', {@myplotx,@myplotfval});

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);

















