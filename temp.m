clear
close all

param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 15;

A = [0,1,0,0,0;
    0,-(param.cf+param.cr)/(param.m*param.vx),(param.cf+param.cr)/(param.m),...
    (-param.lf*param.cf+param.lr*param.cr)/(param.m*param.vx), param.cf/param.m;
    0,0,0,1,0;
    0,(-param.lf*param.cf+param.lr*param.cr)/(param.Iz*param.vx),...
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
dt = 0.05;
sysd = c2d(sys,dt);

Q = [1,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0];

R = [1,0;
    0,0];
N_pred = 100;
mpc = MpcObj(sysd.a,sysd.b,sysd.c,N_pred,Q,R);

x0 = [2,0,0,0,0]';
Yr = zeros(5,N_pred);
mpc.solve(x0,Yr);
path = zeros(N_pred,N_pred*size(B,2));
for k = 1:N_pred
    path(k,k*2) = 1;
end
tic
[U, ~] = mpc.solveConstrained(x0,Yr,[],[],path,zeros(N_pred,1),[-3*ones(1,N_pred);-inf(1,N_pred)],[3*ones(1,N_pred);inf(1,N_pred)]);
toc
U

x0 = [2,0,0,0,0]';
x = zeros(length(x0),N_pred+1);
x(:,1) = x0;
for k = 2:N_pred+1
    x(:,k) = sysd.a*x(:,k-1) + sysd.b*U(:,k-1);
end

plot([x(1,:);x(3,:);x(5,:)]')
% hold on
legend('e1','e2','delta');

t = 0:dt:10;
x0 = [0,2,0,0,0]';
x = zeros(length(x0),length(t));
x(:,1) = x0;
param.delta = 0;
delta_vec = [param.delta];
for k = 2:length(t)
    [U,solve_time(k-1)] = mpc.solveConstrained([x(2:end,k-1);param.delta],Yr,[],[],path,zeros(N_pred,1),[-3*ones(1,N_pred);-inf(1,N_pred)],[3*ones(1,N_pred);inf(1,N_pred)]);
    deltadot = U(1,1);
    param.delta = param.delta + deltadot*dt;
    delta_vec(k) = param.delta;
    zdot = rhs(t(k),x(:,k-1),param);
    x(:,k) = x(:,k-1) + dt*zdot;
end
figure;
plot([x(2,:);x(4,:);delta_vec]')
legend('y','theta','delta');
figure;
plot(x(1,:),x(2,:));
figure;
plot(solve_time);