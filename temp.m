param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 12;

A = [0,1,0,0,0;
    0,-(param.cf+param.cr)/(param.m*param.vx),(param.cf+param.cr)/(param.m),...
    (-param.lf*param.cf+param.lr*param.cr)/(param.m*param.vx), param.cf/param.m;
    0,0,0,1,0;
    0,-(param.lf*param.cf+param.lr*param.cr)/(param.Iz*param.vx),...
    (param.lf*param.cf-param.lr*param.cr)/(param.Iz),...
    -(param.lf^2*param.cf+param.lr^2*param.cr)/(param.Iz*param.vx),param.lf*param.cf/param.Iz;
    0,0,0,0,0];

B = [0,0;
    0,-(param.lf*param.cf-param.lr*param.cr)/param.m - param.vx^2;
    0,0;
    0, -(param.lf^2*param.cf+param.lr^2*param.cr)/(param.Iz);
    1,0];

C = eye(5);

sys = ss(A,B,C,[]);
sysd = c2d(sys,0.05);

Q = [1,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0];

R = [0,0;
    0,0];
N_pred = 50;
mpc = MpcObj(sysd.a,sysd.b,sysd.c,N_pred,Q,R);

x0 = [1,0,0.1,0,0]';
Yr = zeros(5,N_pred);
mpc.solve(x0,Yr);
path = zeros(N_pred,N_pred*size(B,2));
for k = 1:N_pred
    path(k,k*2) = 1;
end
U = mpc.solveConstrained(x0,Yr,[],[],path,zeros(N_pred,1),[-1*ones(1,N_pred);-inf(1,N_pred)],[1*ones(1,N_pred);inf(1,N_pred)]);
U
x = zeros(length(x0),N_pred+1);
x(:,1) = x0;
for k = 2:N_pred+1
    x(:,k) = sysd.a*x(:,k-1) ;%+ sysd.b*U(:,k-1);
end

plot(x')
legend('e1','e1d','e2','e2d','delta');


    