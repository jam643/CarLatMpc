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