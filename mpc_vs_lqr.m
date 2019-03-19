function mpc_vs_lqr()
y = []; x = [];
dt = 0.05;
m = 10;
N_pred = 80;

A = [1, dt; 0, 1];
B = [dt^2/m/2,dt^2/m/2;dt/m,dt/m];
C = [1,0;0,1];

p = size(B,2);
    
Q = [1, 0;
     0, 0];
R = [1e-4, 0;
    0, 2e-4];

x0 = [0; 0];

Xr = ones(1,N_pred);
% Xr = linspace(0,5,N_pred);
Vr = zeros(1,N_pred);
Yr = [Xr;Vr];

[~,~,K] = dare(A,B,Q,R);

mpc_obj = MpcObj(A,B,C,N_pred,Q,R);

U = mpc_obj.solve(x0,Yr);
% U = mpc_obj.solveConstrained(x0,Yr,[],[],[],[],-30*ones(N_pred),30*ones(N_pred));


% plot
t = linspace(0,(N_pred)*dt,N_pred+1);
x(1:2,1) = x0;
y(1:2,1) = C*x0;
x_lqr(1:2,1) = x0;
u_lqr = zeros(p,N_pred);
for k = 2:N_pred+1
    x(:,k) = A*x(:,k-1) + B * U(:,k-1);
    y(:,k) = C*x(:,k);

    u_lqr(:,k-1) =  -K * (x_lqr(:,k-1)-Yr(:,k-1));
    x_lqr(:,k) = A*x_lqr(:,k-1) + B * u_lqr(:,k-1);
end
figure();
set(gcf,'units','normalized','position',[0 0 0.5 1]);
subplot(2,1,1);
plot(t,x(1,:),'k.')
hold on
plot(t,x(2,:),'r.')
plot(t,x_lqr(1,:),'k^')
plot(t,x_lqr(2,:),'r^')
plot(t(2:end),Xr,'k')
plot(t(2:end),Vr,'r')
title('Tracking');
lgnd= legend('x mpc','v mpc','x lqr','v lqr','x setpnt', 'v setpnt');
set(lgnd,'location','best');
subplot(2,1,2);
for k = 1:p
    l1 = plot(t(1:end-1),U(k,:),'k.');
    hold on
    l2 = plot(t(1:end-1),u_lqr(k,:),'r^');
end
lgnd = legend([l1,l2],'mpc','lqr');
set(lgnd,'location','best');
title('Control Effort');
xlabel('Time [s]')
ylabel('Force [N]')
end

