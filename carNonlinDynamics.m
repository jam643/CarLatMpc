function carNonlinDynamics()
t_end = 10;
fps = 100;
N = t_end*fps;
t_array = linspace(0, t_end, N);

param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 5;
param.delta = 0.2;

x0 = 0;
y0 = 0;
vy0 = 0;
theta0 = 0;
thetadot0 = 0;
z0 = [x0,y0,vy0,theta0,thetadot0]';
z_array = zeros(length(z0),N);

f = figure('units','normalized','outerposition',[0 0 1 1],'color',[0.8,0.8,0.8]);
hold on;
set(gca, 'Projection','perspective')
axis equal;
axis vis3d;
set(gca, 'color',[0.85,0.85,0.85],'ZTick', [],'yticklabel',[],'xticklabel',[],'zticklabel',[],'ticklength',[0,0]);
grid on;
% set(gcf,'Renderer','opengl');
axis([0,20,0,20]);
car_handle = fill(0,0,[0.5,0.5,0.5]);
w = 1;
chassis = [-param.lr, param.lf, param.lf, -param.lr;
            w, w, -w, -w];
window = 15;
tic
for k = 1:N-1
    % stop simulation if figure has been closed
    if ~ishandle(f)
        break;
    end
    z_array(:,k+1) = z_array(:,k) + rhs(t_array(k),z_array(:,k),param)*toc;
    tic
    theta = z_array(4,k+1);
    Rot = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    chassisR = Rot*chassis;
    axis([z_array(1,k+1)-window,z_array(1,k+1)+window, z_array(2,k+1)-window, z_array(2,k+1)+window,0,1e-4]);
    view([z_array(4,k+1)*180/pi-90,30]);
    set(car_handle,'xdata',z_array(1,k+1)+chassisR(1,:),'ydata',z_array(2,k+1)+chassisR(2,:))
    pause(0.01);
end
end



function [zdot] = rhs(t,z,param)
x = z(1); y = z(2); vy = z(3);
theta = z(4); thetadot = z(5);

term1 = -param.cf*(atan((vy+param.lf*thetadot)/param.vx) - param.delta)*cos(param.delta);
term2 = param.cr*atan((vy-param.lr*thetadot)/param.vx);
vydot = (term1 - term2)/param.m - param.vx*thetadot;
thetadotdot = (param.lf*term1 + param.lr*term2)/param.Iz;

xdot = param.vx*cos(theta) + vy*sin(theta);
ydot = - vy*cos(theta) + param.vx*sin(theta);

zdot = [xdot, ydot, vydot, thetadot, thetadotdot]';
end