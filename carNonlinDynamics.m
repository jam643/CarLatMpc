function carNonlinDynamics()
t_end = 200;
fps = 100;
N = t_end*fps;
t_array = linspace(0, t_end, N);

param.m = 1915;  %Mass [kg]
param.Iz = 4235;  %Yaw inertia [kg*m^2]
param.lf = 1.453;  %Distance from CG to front axle [m]
param.lr = 1.522;  %Distance from CG to rear axle [m]
param.cf = 90000;  %Front cornering stiffness [N/rad]
param.cr = 116000;  %Rear cornering stiffness [N/rad]
param.vx = 12;
param.delta = 0;
param.max_steer = 1;

leftarrow = false;
rightarrow = false;

x0 = 0;
y0 = 0;
vy0 = 0;
theta0 = 0;
thetadot0 = 0;
z0 = [x0,y0,vy0,theta0,thetadot0]';
z_array = zeros(length(z0),N);

f = figure('units','normalized','outerposition',[0 0 1 1],'color',[0.8,0.8,0.8]);
b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',0, 'min',-param.max_steer, 'max',param.max_steer);
addlistener(b, 'Value', 'PostSet', @updateSteer);
% b.Callback = @updateSteer; 
% set(f,'KeyPressFcn',@KeyPressCb, 'KeyReleaseFcn',@KeyReleaseCb) ;
hold on;
set(gca, 'Projection','perspective')
axis equal;
axis vis3d;
axis manual;
set(gca, 'color',[0.85,0.85,0.85],'ZTick', [],'yticklabel',[],'xticklabel',[],'zticklabel',[],'ticklength',[0,0]);
grid on;
axis([0,20,0,20,0,1e-6]);
track_handle = plot(0,0,'k');
car_handle = fill(0,0,[0.5,0.5,0.5]);
tire_fl_handle = fill(0,0,[0.2,0.2,0.2]);
tire_fr_handle = fill(0,0,[0.2,0.2,0.2]);
tire_rl_handle = fill(0,0,[0.2,0.2,0.2]);
tire_rr_handle = fill(0,0,[0.2,0.2,0.2]);
steer_handle = fill3(0,0,0,[0.2,0.2,0.2]);
car_width = 1;
tire_width = car_width/8;
tire_len = (param.lr+param.lf)/8;
chassis = [-param.lr-tire_len, param.lf+tire_len, param.lf+tire_len, -param.lr-tire_len;
            car_width+tire_width, car_width+tire_width, -car_width-tire_width, -car_width-tire_width];
tire = [-tire_len, tire_len, tire_len, -tire_len;
            tire_width, tire_width, -tire_width, -tire_width];
steer_len = tire_width;
car_height = car_width/2;
steer = [0,0,0,0;
    -steer_len, steer_len, steer_len, -steer_len;
    steer_len+car_height, steer_len+car_height, -steer_len+car_height, -steer_len+car_height];
window = 15;
cameraAngle = theta0*180/pi;
cameraAngleRate = 0;
cameraAngleErrorCum = 0;
tLastUpdate = tic;
tStart = tic;
k = 1;
steer_rate = 3;
steer_return_rate = 3;
max_steer = 0.6;
while ishandle(f) 
%     param.delta = 0.5*sin(toc(tStart)*2*pi/5);
    z_array(:,k+1) = z_array(:,k) + rhs(t_array(k),z_array(:,k),param)*toc(tLastUpdate);
    cameraAngleError = z_array(4,k+1)*180/pi - cameraAngle;
    cameraAngleRate = cameraAngleRate + toc(tLastUpdate)*(10*cameraAngleError-6*cameraAngleRate);
    cameraAngle = cameraAngle + toc(tLastUpdate)*cameraAngleRate;
%     cameraAngle = 180*z_array(4,k+1)/pi;
    
%     if rightarrow || leftarrow
%         if leftarrow
%             param.delta = min(max_steer,param.delta + steer_rate*toc(tLastUpdate));
%         end
%         if rightarrow
%             param.delta = max(-max_steer,param.delta - steer_rate*toc(tLastUpdate));
%         end
%     else
%         param.delta = param.delta - sign(param.delta)*steer_return_rate*toc(tLastUpdate);
%     end
    
    tLastUpdate = tic;
    theta = z_array(4,k+1);
    
    set(track_handle,'xdata',z_array(1,1:k+1),'ydata',z_array(2,1:k+1));
    chassisR = homogTransform(theta,z_array(1,k+1),z_array(2,k+1),chassis);
    tire_fl_rot = homogTransform(theta,z_array(1,k+1),z_array(2,k+1),homogTransform(param.delta,param.lf,car_width,tire));
    tire_fr_rot = homogTransform(theta,z_array(1,k+1),z_array(2,k+1),homogTransform(param.delta,param.lf,-(car_width),tire));
    tire_rl_rot = homogTransform(theta,z_array(1,k+1),z_array(2,k+1),homogTransform(0,-param.lr,car_width,tire));
    tire_rr_rot = homogTransform(theta,z_array(1,k+1),z_array(2,k+1),homogTransform(0,-param.lr,-(car_width),tire));
    axis([z_array(1,k+1)-window,z_array(1,k+1)+window, z_array(2,k+1)-window, z_array(2,k+1)+window]);
    view([cameraAngle-90,30]);
    set(car_handle,'xdata',chassisR(1,:),'ydata',chassisR(2,:))
    set(tire_fl_handle,'xdata',tire_fl_rot(1,:),'ydata',tire_fl_rot(2,:))
    set(tire_fr_handle,'xdata',tire_fr_rot(1,:),'ydata',tire_fr_rot(2,:))
    set(tire_rl_handle,'xdata',tire_rl_rot(1,:),'ydata',tire_rl_rot(2,:))
    set(tire_rr_handle,'xdata',tire_rr_rot(1,:),'ydata',tire_rr_rot(2,:))
   
    pause(0.02)
    k = k + 1;
end

function KeyPressCb(~,evnt)
if strcmpi(evnt.Key,'leftarrow')
    leftarrow = true;
elseif strcmpi(evnt.Key,'rightarrow')
    rightarrow = true;
end
end

function KeyReleaseCb(~,evnt)
if strcmpi(evnt.Key,'leftarrow')
    leftarrow = false;
elseif strcmpi(evnt.Key,'rightarrow')
    rightarrow = false;
end
end

function updateSteer(src,event)
    param.delta = -event.AffectedObject.Value;
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

function [shapeRot] = homogTransform(theta,x,y,shape)
Rot = [cos(theta),-sin(theta);sin(theta),cos(theta)];
shapeRot = Rot*shape + repmat([x;y],1,length(shape));
end