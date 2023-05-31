clear variables;
%% initializing for parameters
% a（m）：质心到前轴的距离
% b（m）：质心到后轴的距离

m = 1270;
mf = 71*2;
mr = 71*2;

g = 9.8;

Izz = 1536.7;
Ixx = 536.6;
Ixz = 0;

a = 1.015;
b = 2.910-a;
hb = 0.115;
lf = 0;
lr = 0;

Crll = -3577*1.575*1.575/2;
Krll = -(27000*1.575*1.575/2);

%% 调参数 
% Kty1（N/rad）：左前轮侧偏刚度
Kty1 =(-4600*15*2);
Kty2 =(-1800*15*2);

%% basic simulation parameters 仿真时间和步长
TimeS = 20;
step = 0.01;
Time = 0:step:TimeS;
nstep = length(Time);

%% initializing for variables
% δ （rad）：前轮转角
% delta = 2*pi/180 * sin(2*pi/10 * Time);%%%%sine input
delta =15 * pi/180+zeros(nstep,1);%%%%step input

u = 10;
v = zeros(nstep,1);v(1) = 0;
p = zeros(nstep,1);p(1) = 0;
r = zeros(nstep,1);r(1) = 0;
dv = zeros(nstep,1);dv(1) = 0;
dp = zeros(nstep,1);dp(1) = 0;
dr = zeros(nstep,1);dr(1) = 0;

%(m/s^2)
rll = zeros(nstep,1);
drll = zeros(nstep,1);drll(1) = 0;

alpha1 = zeros(nstep,1);
alpha2 = zeros(nstep,1);

fy1 = zeros(nstep,1);
fy2 = zeros(nstep,1);

%% adding global variables
dyaw = zeros(nstep,1);
dpth = zeros(nstep,1);
yaw = zeros(nstep,1);

x = zeros(nstep,1);
y = zeros(nstep,1);
z = zeros(nstep,1);
dx = zeros(nstep,1);
dy = zeros(nstep,1);
dz = zeros(nstep,1);
%% integration 
for i=1:(nstep-1)
    
    [dyaw(i),trash1,drll(i)]=AngleRateA(p(i),0,r(i),rll(i),0);
    %-------------算法需改进-----------------------------------------
    rll(i+1)=rll(i)+step*drll(i);
    yaw(i+1)= yaw(i)+step*dyaw(i);

    %1大地-车身坐标变换矩阵
    Cga = CgaCal(rll(i),0,yaw(i));
    temp1 = Cga*[u;v(i);0];
    dx(i) = temp1(1);
    dy(i) = temp1(2);
    dz(i) = temp1(3);
    %------------算法需改进------------------------------------
    x(i+1) = x(i) + step * dx(i);
    y(i+1) = y(i) + step * dy(i);
    z(i+1) = z(i) + step * dz(i);
    
    beta = atan(v(i)/u);
    alpha1(i) = beta + a*r(i)/u - delta(i);
    alpha2(i) = beta - b*r(i)/u;
    
    if(abs(u)<2.5 && u~=0)   %低速下限制侧偏角的值以提高稳定性
        alpha1(i) = alpha1(i) * ((1-0.01)*abs(u)/2.5 + 0.01);
        alpha2(i) = alpha2(i) * ((1-0.01)*abs(u)/2.5 + 0.01);
    end
    
    fy1(i)=Kty1*alpha1(i);
    fy2(i)=Kty2*alpha2(i);
    
    M = [1,0,0,0;
        0,(m + mf + mr),m * hb,(a*mf - b*mr);
        0,(a*mf - b*mr),Ixz,Izz;
        0,m * hb,Ixx,Ixz];
    
    FF = [p(i);
          (fy1(i)*cos(delta(i)) + fy2(i))*sign(u) - (m + mf + mr) * u * r(i);
          (a*fy1(i)*cos(delta(i)) - b*fy2(i))*sign(u) - (a*mf - b*mr) * u * r(i);
          lf * fy1(i) + lr * fy2(i) + Crll * p(i) + (Krll - m * g * hb) * rll(i)- m * hb *u*r(i);];
      
    
      
    temp = M\FF;
    
    drll(i) = temp(1);
    dv(i) = temp(2);
    dp(i) = temp(3);
    dr(i) = temp(4);
    
    rll(i+1) = rll(i) + step * drll(i);
    v(i+1) = v(i) + step * dv(i);
    p(i+1) = p(i) + step * dp(i);
    r(i+1) = r(i) + step * dr(i);
    
end

%% post processing
figure(10);
subplot(4,4,1);plot(Time,v);title('v（m/s）');
subplot(4,4,2);plot(Time,57.3*r);title('r(deg/s)');
subplot(4,4,3);plot(Time,dv);title('dv(m/s2)');
subplot(4,4,4);plot(Time,57.3*dr);title('dr(deg/s2)');
subplot(4,4,5);plot(Time,57.3*alpha1);title('alpha1(deg)');
subplot(4,4,6);plot(Time,57.3*alpha2);title('alpha2(deg)');
subplot(4,4,7);plot(Time,fy1);title('fyF(N)');
subplot(4,4,8);plot(Time,fy2);title('fyR(N)');
subplot(4,4,9);plot(Time,57.3*p);title('p(deg/s)');
subplot(4,4,10);plot(Time,57.3*dp);title('dp(deg/s2)');
subplot(4,4,11);plot(Time,57.3*rll);title('rll(deg)');
subplot(4,4,12);plot(Time,57.3*drll);title('drll(deg/s)');
subplot(4,4,13);plot(Time,x);title('x(m)');
subplot(4,4,14);plot(Time,y);title('y(m)');
subplot(4,4,15);plot(Time,z);title('z(m)');
subplot(4,4,16);plot(x,y);title('轨迹(m)');