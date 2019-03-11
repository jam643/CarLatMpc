param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 12;

A = [0,1,0,0,0;
    0,-(param.cf+param.cr)/(param.m*param.vx),(param.cf+param.cr)/(param.m),
    (-param.lf*param.cf+parm.lr*param.cr)/(param.m*param.vx), param.cf/param.m;
    0,0,0,1,0;
    0,-(param.lf*param.cf+parm.lr*param.cr)/(param.Iz*param.vx),
    (param.lf*param.cf-parm.lr*param.cr)/(param.Iz),
    -(param.lf^2*param.cf+parm.lr^2*param.cr)/(param.Iz*param.vx),parm.lf*param.cf/param.Iz;
    0,0,0,0,0];

B = [0,0;
    0,-(param.lf*param.cf-param.lr*param.cr)/param.m - param.vx^2;
    0,0;
    0, -(param.lf^2*param.cf+parm.lr^2*param.cr)/(param.Iz);
    1,0];
    