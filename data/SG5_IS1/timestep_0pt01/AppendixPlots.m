



%--------------------------------------------------------------------------
% Calculates and plots u, v, and vel velocity time histories, and
% calculates u', v', and vel' accelerations using the central difference
% method.

% Loads x and y position, and u and v velocity, and time
open('u_vs_time_fit.mat');
u_vs_time_fit=ans.u_vs_time_fit;
u=u_vs_time_fit.yfit;
time=u_vs_time_fit.xi;
[J,I]=size(u);
%u=u/1000;       %to convert from mm/s to m/s

open('v_vs_time_fit.mat');
v_vs_time_fit=ans.v_vs_time_fit;
v=v_vs_time_fit.yfit;
%v=v/1000;       %to convert from mm/s to m/s

% generate vel array (=sqrt(v^2)+(v^2))
vel=(sqrt((u.^2)+(v.^2)));

% use central difference method to derive acceleration data from velocities
t=2;            % start at row 1 of WL array

for j=1:(J-2)
    u_prime(j)=(u(t+1)-u(t-1))/(time(t+1)-time(t-1));
    v_prime(j)=(v(t+1)-v(t-1))/(time(t+1)-time(t-1));
    vel_prime(j)=(vel(t+1)-vel(t-1))/(time(t+1)-time(t-1));
t=t+1;
end

u_max=max(u);
v_max=min(v);
vel_max=max(vel);
u_v_vel_max=[u_max v_max vel_max];

u=u./1000;
v=v./1000;
vel=vel./1000;
u_prime=u_prime./1000;
v_prime=v_prime./1000;
vel_prime=vel_prime./1000;

figure(1)
plot(time,vel,'k:o');
title('Landslide Velocity Time History','fontsize',12);
xlabel('time (s)');
ylabel('velocity (m/s)');
axis([0 6 0 1.8]);
grid on;

figure(2)
plot(time(2:(J-1)),vel_prime,'k:o');
title('Landslide Acceleration Time History','fontsize',12);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
axis([0 6 -2.500 2]);
grid on;


%--------------------------------------------------------------------------
open('PosAll_1_31.mat');
Zall_Xmod=ans.Zall_Xmod;
Zall_WL=ans.Zall_WL/1000;
Zall_T=ans.Zall_T;

y1=Zall_WL(25,(120:end));
y2=Zall_WL(40,(120:end));
y3=Zall_WL(55,(120:end));
y4=Zall_WL(70,(120:end));
y5=Zall_WL(85,(120:end));
y6=Zall_WL(100,(120:end));

Zall_Xmod=Zall_Xmod./1000;

figure(3)
set(axes('FontSize', 12));
subplot(6,1,1);plot(Zall_Xmod(120:end),y1,'k');axis([0 10 -0.04 0.04]);grid on;
ylabel('eta (m)');title('Water Surface Profile (time =0.600sec)');
set(gca,'XTickLabel',{''})

subplot(6,1,2);plot(Zall_Xmod(120:end),y2,'k');axis([0 10 -0.04 0.04]);grid on;
ylabel('eta (m)');title('(time =1.600sec)')
set(gca,'XTickLabel',{''})

subplot(6,1,3);plot(Zall_Xmod(120:end),y3,'k');axis([0 10 -0.04 0.04]);grid on;
ylabel('eta (m)');title('(time =2.600sec)')
set(gca,'XTickLabel',{''})

subplot(6,1,4);plot(Zall_Xmod(120:end),y4,'k');axis([0 10 -0.04 0.04]);grid on;
ylabel('eta (m)');title('(time =3.600sec)')
set(gca,'XTickLabel',{''})

subplot(6,1,5);plot(Zall_Xmod(120:end),y5,'k');axis([0 10 -0.04 0.04]);grid on;
ylabel('eta (m)');title('(time =4.600sec)')
set(gca,'XTickLabel',{''})

subplot(6,1,6);plot(Zall_Xmod(120:end),y6,'k');axis([0 10 -0.04 0.04]);grid on;
xlabel('x (m)'); ylabel('eta (m)')
title('(time =5.600sec)')



y1=Zall_WL((16:end),1373);
y2=Zall_WL((16:end),3879);
y3=Zall_WL((16:end),6386);
y4=Zall_WL((16:end),8892);
y5=Zall_WL((16:end),11398);
y6=Zall_WL((16:end),13904);

figure(4)
subplot(6,1,1);plot(Zall_T(16:end),y1,'k');axis([0 7 -0.04 0.04]);grid on;
set(gca,'XTickLabel',{''}); ylabel('eta (m)');title('Water Surface Time History (x = 0.500m)')

subplot(6,1,2);plot(Zall_T(16:end),y2,'k');axis([0 7 -0.04 0.04]);grid on;
set(gca,'XTickLabel',{''}); ylabel('eta (m)');title('(x = 1.500m)')

subplot(6,1,3);plot(Zall_T(16:end),y3,'k');axis([0 7 -0.04 0.04]);grid on;
set(gca,'XTickLabel',{''}); ylabel('eta (m)');title('(x = 2.500m)')

subplot(6,1,4);plot(Zall_T(16:end),y4,'k');axis([0 7 -0.04 0.04]);grid on;
set(gca,'XTickLabel',{''}); ylabel('eta (m)');title('(x = 3.500m)')

subplot(6,1,5);plot(Zall_T(16:end),y5,'k');axis([0 7 -0.04 0.04]);grid on;
set(gca,'XTickLabel',{''}); ylabel('eta (m)');title('(x = 4.500m)')

subplot(6,1,6);plot(Zall_T(16:end),y6,'k');axis([0 7 -0.04 0.04]);grid on;
xlabel('time (s)'); ylabel('eta (m)');title('(x = 5.500m)')
set(gca,'XTickLabel',{''});

%--------------------------------------------------------------------------

% Calculation of the wave potential energy (Ep), and plot against
% time.  

clear

open('PosAll_1_31.mat');
Zall_Xmod=ans.Zall_Xmod;
Zall_WL=ans.Zall_WL;
Zall_T=ans.Zall_T;
Zall_Xmod=0.001*(Zall_Xmod);
Zall_WL=0.001*(Zall_WL);
[J,I]=size(Zall_WL);

roh=1000;           % kg per cubic metre, density of fluid
g=9.81;             % m per sec^2, acceleration of gravity
width=0.25;         % width of flume

% Central Difference integration
x=2;            % start at col 2 of WL array
t=1;            % start at row 1 of WL array

for j=1:J
    for i=1:I-2
        wave_Ep_Array(t,x)=0.5*roh*g*width*((((Zall_WL(t,x+1)+Zall_WL(t,x-1))/2)^2)*((Zall_Xmod(x+1,1)-Zall_Xmod(x-1,1))/2));
    x=x+1;
end;
t=t+1;
x=2;
end
wave_Ep=sum(wave_Ep_Array');

run 'Block_Ek_Ep_vs_T';

figure(5)
plot(Zall_T,wave_Ep,'xk')
title('Wave Potential Energy - 0.5 \int \rho g \eta^2 width dx');
xlabel('time (s)');
ylabel('energy (joules)');
axis([0 8 0 1.5]);   
grid on;

figure(6)
plot(Zall_T,wave_Ep,'k+:')
hold on;
plot(time,Ep,'ko:');
plot(time,Ek,'k*:');
hold off;
title('Wave and Slider Energy' );
xlabel('time (s)');
ylabel('energy (joules)');
axis([0 8 0 30]);  
legend('Wave Ep','Slider Ep','Slider Ek',1);
grid on;


%--------------------------------------------------------------------------

% Calculate positive and negative water level (vs space) envelopes and
% position and time maximum and minimum WL occur

% Loads Position, Time, and Water Level variables from saved combined .mat
% file
open('PosAll_1_31.mat');
Zall_Xmod=ans.Zall_Xmod;
Zall_WL=ans.Zall_WL/1000;
Zall_T=ans.Zall_T;
Zall_Xmod=0.001*(Zall_Xmod);

% define WL array size
[m,n]=size(Zall_WL);

% calculate WL maxima and minima at each spatial location
for i=1:n;
    positiveWL(i)=max(Zall_WL(:,i));
    negativeWL(i)=min(Zall_WL(:,i));
end

% find max and min WL at each time step, and finds its corresponding
% location
Zall_WL_abs=abs(Zall_WL);
for i=1:m;
    [B,maxXmodRow]=max(Zall_WL(i,:));
    [C,minXmodRow]=min(Zall_WL(i,:));
    max_positiveWL(i)=B;
    max_negativeWL(i)=C;
    maxpos(i)=Zall_Xmod(maxXmodRow);
    minpos(i)=Zall_Xmod(minXmodRow);
    max_amplitude(i)=max(Zall_WL_abs(i,:));
end

figure(7)
plot(Zall_Xmod,positiveWL,'k');
hold on;
plot(Zall_Xmod,negativeWL,'k');
hold off;
title('Maximum and Minimum Water Level Envelope');
xlabel('x (m)');
ylabel('eta (m)');
axis([0 10.200 -0.04 0.04]);
grid on;

figure(8)
plot(Zall_T,max_positiveWL,'k');
hold on;
plot(Zall_T,max_negativeWL,'k');
hold off;
title('Maximum and Minimum Water Level Time History');
xlabel('time (s)');
ylabel('eta (m)');
axis([0 8 -0.04 0.04]);
grid on;

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------



%--------------------------------------------------------------------------