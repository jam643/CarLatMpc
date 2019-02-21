function forcesProTest()

%% NLP relevant dimensions
nlp.N = 50; % horizon length
nlp.nvar = 6; % number of stage variables
nlp.neq = 5; % number of equality constraints
nlp.nh = 0; % number of inequality constraints
nlp.npar = 0; % number of runtime parameters

%% Objective function
nlp.objective = @objective;

%% Matrix equality constraints
nlp.eq = @dynamics;
nlp.E = [zeros(5,1),eye(5)]; %why?

%% Inequality constraints
nlp.lb = [-inf,-inf,-inf,-inf,-inf,-1];
nlp.ub = [inf,inf,inf,inf,inf,1];

%% Init condition
nlp.xinitidx = 1:5;


% Get the default solver options
codeoptions = getOptions('FORCESNLPsolver');
% Generate solver
status = FORCES_NLP (nlp , codeoptions );
end

function f = objective(w)
x = w(1);
u = w(6);
f = x^2 + u^2;
end

function [zdot] = continuous_dynamics(z,u)
%% Parameters
param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 5;
delta = u(1);

x = z(1); y = z(2); vy = z(3);
theta = z(4); thetadot = z(5);

term1 = -param.cf*(atan((vy+param.lf*thetadot)/param.vx) - delta)*cos(delta);
term2 = param.cr*atan((vy-param.lr*thetadot)/param.vx);
vydot = (term1 - term2)/param.m - param.vx*thetadot;
thetadotdot = (param.lf*term1 + param.lr*term2)/param.Iz;

xdot = param.vx*cos(theta) + vy*sin(theta);
ydot = - vy*cos(theta) + param.vx*sin(theta);

zdot = [xdot, ydot, vydot, thetadot, thetadotdot]';
end

function znext = dynamics(w)
z = zeros(5,1);
z = w(1:5);
u = zeros(1,1);
u = w(6);
h = 0.05; %step size
znext = RK4(z,u,@continuous_dynamics,h);
end