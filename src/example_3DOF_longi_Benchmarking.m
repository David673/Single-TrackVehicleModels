clear variables;
%% initializing for parameters
% a(m):The distance from the center of mass to the front axle
% b(m):The distance from the center of mass to the rear axle

a = 1.015;
b = 2.910-a;
Izz = 1536.7;
m = 1270;

%% Adjust the parameters 
% Kty1(N/rad):Side-slip stiffness of the left front wheel
Kty1 =(-4600*15*2);
Kty2 =(-1800*15*2);

%% basic simulation parameters
TimeS = 80;
step = 0.0005;
Time = 0:step:TimeS;
nstep = length(Time);

%% initializing for variables
% delta(rad):steering angle of the front wheels

delta = 15*pi/180 * sin(2*pi/10 * Time);%%%%sine input
% delta =15 * pi/180+zeros(nstep,1);%%%%step input

u = zeros(nstep,1);u(1) = 10;
v = zeros(nstep,1);v(1) = 0;
r = zeros(nstep,1);r(1) = 0;

%(m/s^2)
du = zeros(nstep,1);du(1) = 0;
dv = zeros(nstep,1);dv(1) = 0;
dr = zeros(nstep,1);dr(1) = 0;

alpha1 = zeros(nstep,1);
alpha2 = zeros(nstep,1);

fx1 = zeros(nstep,1);
fx2 = zeros(nstep,1);
fy1 = zeros(nstep,1);
fy2 = zeros(nstep,1);

% %% additional variables
% yaw = zeros(nstep,1);
% x = zeros(nstep,1);
% y = zeros(nstep,1);
% z = zeros(nstep,1);
% dx = zeros(nstep,1);
% dy = zeros(nstep,1);
% dz = zeros(nstep,1);

%% integration 
for i=1:(nstep-1)
%     yaw(i+1)= yaw(i)+step*r(i);
%     %
%     Cga = CgaCal(0,0,yaw(i));
%     temp1 = Cga*[u(i);v(i);0];
%     dx(i) = temp1(1);
%     dy(i) = temp1(2);
%     dz(i) = temp1(3);
%     
%     x(i+1) = x(i) + step * dx(i);
%     y(i+1) = y(i) + step * dy(i);
%     z(i+1) = z(i) + step * dz(i);
 
    
    beta = atan(v(i)/u(i));
    alpha1(i) = beta + a*r(i)/u(i) - delta(i);
    alpha2(i) = beta - b*r(i)/u(i);
    
%     if(abs(u(i))<2.5 && u(i)~=0)   
%         alpha1(i) = alpha1(i) * ((1-0.01)*abs(u(i))/2.5 + 0.01);
%         alpha2(i) = alpha2(i) * ((1-0.01)*abs(u(i))/2.5 + 0.01);
%     end
    
    fy1(i)=Kty1*alpha1(i);
    fy2(i)=Kty2*alpha2(i);
    
    dr(i) = 1/Izz * (a*fy1(i) * cos(delta(i)) - b*fy2(i));
    dv(i) = -u(i) * r(i) + 1/m * (fy1(i) * cos(delta(i)) + fy2(i));
    du(i) = r(i) * v(i) + 1/m * (fx1(i)+fx2(i) - fy1(i)*sin(delta(i)));
    
%     dr(i) = 1/Izz * (a*fy1(i) - b*fy2(i));
%     dv(i) = -u(i) * r(i) + 1/m * (fy1(i) + fy2(i));
%     du(i) = r(i) * v(i) + 1/m * (fx1(i) + fx2(i));
    
    r(i+1) = r(i) + step * dr(i);
    v(i+1) = v(i) + step * dv(i);
    u(i+1) = u(i) + step * du(i);
end

% titleS={'Time','delta','alpha1','alpha2','fx1','fx2','fy1','fy2','u','v','r','du','dv','dr'};
% dataOut = table(Time',delta',alpha1,alpha2,fx1,fx2,fy1,fy2,u,v,r,du,dv,dr,'VariableNames',titleS);
% writetable(dataOut,'dataOut3DOF.csv');
% 
% figure(13);
% subplot(4,4,1);plot(Time,v);title('v（m/s）');
% subplot(4,4,2);plot(Time,57.3*r);title('r(deg/s)');
% subplot(4,4,3);plot(Time,dv);title('dv(m/s2)');
% subplot(4,4,4);plot(Time,57.3*dr);title('dr(deg/s2)');
% subplot(4,4,5);plot(Time,57.3*alpha1);title('alpha1(deg)');
% subplot(4,4,6);plot(Time,57.3*alpha2);title('alpha2(deg)');
% subplot(4,4,7);plot(Time,fy1);title('fyF(N)');
% subplot(4,4,8);plot(Time,fy2);title('fyR(N)');
% subplot(4,4,9);plot(Time,u);title('u（m/s）');
% subplot(4,4,10);plot(Time,du);title('du（m/s2）');
% subplot(4,4,11);plot(Time,fx1);title('fxF（N）');
% subplot(4,4,12);plot(Time,fx2);title('fxR（N）');
% % subplot(4,4,13);plot(Time,x);title('x（m）');
% % % subplot(4,4,14);plot(Time,y);title('y（m）');
% % subplot(4,4,15);plot(Time,z);title('z（m）');
% % subplot(4,4,16);plot(x,y);title('Traj(m');
% 
% fy1sin = zeros(1,nstep);
% for i=1:nstep
%     fy1sin(i)=fy1(i)*sin(delta(i));
% end
% figure(14);
% plot(Time,fy1sin);title('fy1*sin(delta)(N)');