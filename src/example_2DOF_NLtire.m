clear variables;
close all;
%% initializing for parameters
% a(m):The distance from the center of mass to the front axle
% b(m):The distance from the center of mass to the rear axle

a = 1.615;
b = 2.910-a;
Izz = 1536.7;
m = 1270;
g = 9.8;

%% Adjust the parameters 
% Kty1(N/rad):Side-slip stiffness of the left front wheel
Kty1 =(-800*15*2);
Kty2 =(-800*15*2);

%% basic simulation parameters 
TimeS = 30;
step = 0.01;
Time = 0:step:TimeS;
nstep = length(Time);

%% initializing for variables
% delta(rad):steering angle of the front wheels

% delta = 2*pi/180 * sin(2*pi/10 * Time);%%%%sine input
delta = 15 * pi/180+zeros(nstep,1);%%%%step input

u = 20;
v = zeros(nstep,1);v(1) = 0;

r = zeros(nstep,1);r(1) = 0;

%(m/s^2)
dv = zeros(nstep,1);dv(1) = 0;
dr = zeros(nstep,1);dr(1) = 0;

alpha1 = zeros(nstep,1);
alpha2 = zeros(nstep,1);

fy1 = zeros(nstep,1);
fy2 = zeros(nstep,1);

%% additional variales
yaw = zeros(nstep,1);

x = zeros(nstep,1);
y = zeros(nstep,1);
z = zeros(nstep,1);
dx = zeros(nstep,1);
dy = zeros(nstep,1);
dz = zeros(nstep,1);



%% integration 
for i=1:(nstep-1)
    
    yaw(i+1)= yaw(i)+step*r(i);
    %Transformation matrix between ground and vehicle body
    Cga = CgaCal(0,0,yaw(i));
    temp1 = Cga*[u;v(i);0];
    dx(i) = temp1(1);
    dy(i) = temp1(2);
    dz(i) = temp1(3);
    
    x(i+1) = x(i) + step * dx(i);
    y(i+1) = y(i) + step * dy(i);
    z(i+1) = z(i) + step * dz(i);

    
    beta = atan(v(i)/u);
    alpha1(i) = beta + a*r(i)/u - delta(i);
    alpha2(i) = beta - b*r(i)/u;
    
    if(abs(u)<2.5 && u~=0)   %Limit the alpha angles at low speeds to improve stability
        alpha1(i) = alpha1(i) * ((1-0.01)*abs(u)/2.5 + 0.01);
        alpha2(i) = alpha2(i) * ((1-0.01)*abs(u)/2.5 + 0.01);
    end
    
    fy1(i)=Kty1*alpha1(i);
    fy2(i)=Kty2*alpha2(i);
    if (fy1(i)>=0.9*m*g*(b/(a+b)))
        fy1(i) = 0.9*m*g*(b/(a+b));
    end
    if (fy2(i)>=0.9*m*g*(a/(a+b)))
        fy2(i) = 0.9*m*g*(a/(a+b));
    end
    
%     dr(i) = 1/Izz * (a*fy1(i) - b*fy2(i));
%     dv(i) = -u * r(i) + 1/m * (fy1(i) + fy2(i));
    
    dr(i) = 1/Izz * (a*fy1(i) * cos(delta(i)) - b*fy2(i))*sign(u);
    dv(i) = -u * r(i) + 1/m * (fy1(i) * cos(delta(i)) + fy2(i))*sign(u);
    
    r(i+1) = r(i) + step * dr(i);
    v(i+1) = v(i) + step * dv(i);
end

figure(8);
subplot(3,4,1);plot(Time,v);title('v£¨m/s£©');
subplot(3,4,2);plot(Time,57.3*r);title('r(deg/s)');
subplot(3,4,3);plot(Time,dv);title('dv(m/s2)');
subplot(3,4,4);plot(Time,57.3*dr);title('dr(deg/s2)');
subplot(3,4,5);plot(Time,57.3*alpha1);title('alpha1(deg)');
subplot(3,4,6);plot(Time,57.3*alpha2);title('alpha2(deg)');
subplot(3,4,7);plot(Time,fy1);title('fyF(N)');
subplot(3,4,8);plot(Time,fy2);title('fyR(N)');
subplot(3,4,9);plot(Time,57.3*yaw);title('yaw(deg)');
subplot(3,4,10);plot(Time,x);title('x(m)');
subplot(3,4,11);plot(Time,y);title('y(m)');
subplot(3,4,12);plot(x,y);title('Trajectory(m)');

figure(9);
subplot(1,2,1);plot(Time,fy1),title('Fy1(N)');
subplot(1,2,2);plot(Time,fy2),title('Fy2(N)');

figure(10);
plot(x,y);hold on;axis equal;
for i=1:nstep
    if mod(i,round(nstep/40))==0
        quiver(x(i),y(i),10*cos(yaw(i)),10*sin(yaw(i)),1,'r','LineWidth',1);
    end
end
