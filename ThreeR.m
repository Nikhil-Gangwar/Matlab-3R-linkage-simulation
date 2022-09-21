clc;
clear all;

m1=1;
m2=1;
m3=1;

l1=1;
l2=0.75;
l3=0.5;

g=9.8;

th_1(1)=pi/3;
thdot_1(1)=0;
th_2(1)=pi/4;
thdot_2(1)=0;
th_3(1)=pi/2;
thdot_3(1)=0;

ts=10;
dt=0.05;
t=0:dt:10;


for i=1:length(t)
    
    sig1=3*cos(th_2(i))/2 + cos(th_3(i)) + cos(th_2(i)+th_3(i))/2 + 5/3;
    sig2=cos(th_3(i))/2 + cos(th_2(i) + th_3(i))/2 +1/3;
    M=[3*cos(th_2(i))+cos(th_3(i))+cos(th_2(i)+th_3(i))+4 sig1 sig2; sig1 cos(th_3(i))+5/3 cos(th_3(i))/2 + 1/3; sig2 cos(th_3(i))/2 + 1/3 1/3];


    sig1H=(981*cos(th_1(i)) + th_2(i) + th_3(i))/200;
    H=[(981*cos(th_1(i)))/40 + sig1H + (2943*cos(th_1(1) + th_2(i)))/200 ; (991*cos(th_1(i) + th_2(i) + th_3(i)))/200 + (2973*cos(th_1(i) + th_2(i)))/200 ; sig1H]; 

   %     
    sig5C=(thdot_1(i)+thdot_2(i)+thdot_3(i))^2;
    sig6C=th_1(i)+th_2(i)+th_3(i);
    sig7C=cos(th_1(i))*(thdot_1(i)^2);
    sig8C=sin(th_1(i))*(thdot_1(i)^2);
    sig9C=(thdot_1(i)+(thdot_2(i))^2);
    sig10C=th_1(i)+th_2(i);
    sig1C=sig6C + sig5C*sin(sig6C)/2 + sig9C*sin(sig10C);
    sig2C= sig5C*cos(sig6C)/2 + sig9C*cos(sig10C) + sig1C;
    sig3C=sig9C*cos(sig10C)/2 + sig1C;
    sig4C=sig8C + (sig9C*sin(sig10C))/2;
%

    C=[sig2C*(sin(th_1(i)) + sin(sig6C)/2 + sin(sig10C)) - (cos(th_1(i)) + cos(sig6C)/2 + cos(sig10C))*sig1C - (cos(th_1(i)) + cos(sig10C)/2)*sig4C + (sin(th_1(i)) + sin(sig10C)/2)*sig3C ; (sin(sig6C)/2 + sin(sig10C))*sig2C - (cos(sig6C)/2 + cos(sig10C))*sig1C + (sin(sig10C)*sig3C)/2 - (cos(sig10C)*sig4C)/2  ; (sin(sig6C)*sig2C)/2 - (cos(sig6C)*sig1C)/2];
    
    Minv=inv(M);
    thddot=-Minv*(C+H);
    thdot_1(i+1)=thdot_1(i) + thddot(1)*dt;
    th_1(i+1)=th_1(i)+thdot_1(i+1)*dt;
    thdot_2(i+1)=thdot_2(i) + thddot(1)*dt;
    th_2(i+1)=th_2(i)+thdot_2(i+1)*dt;
    thdot_3(i+1)=thdot_3(i) + thddot(2)*dt;
    th_3(i+1)=th_3(i)+thdot_3(i+1)*dt;

    hold on;
    plot([0 l1*cos(th_1(i))],[0 l1*sin(th_1(i))],'b-o');
    plot([l1*cos(th_1(i)) l1*cos(th_1(i))+l2*cos(th_2(i))],[l1*sin(th_1(i)) l1*sin(th_1(i))+l2*sin(th_2(i))],'r-o');
    plot([l1*cos(th_1(i))+l2*cos(th_2(i)) l1*cos(th_1(i))+l2*cos(th_2(i))+l3*cos(th_3(i))],[l1*sin(th_1(i))+l2*sin(th_2(i)) l1*sin(th_1(i))+l2*sin(th_2(i))+l3*sin(th_3(i))],'g-o');
    axis([-2 2 -2 2])
    pause(dt+0.1);
    hold off;


end