% Calculates and plots u, v, and vel velocity time histories, and
% calculates u', v', and vel' accelerations using the central difference
% method.

% Loads x and y position, and u and v velocity, and time
open('u_vs_time_fit.mat');
u_vs_time_fit=ans.u_vs_time_fit;
u=u_vs_time_fit.yfit;
time=u_vs_time_fit.xi
[J,I]=size(u);
%u=u/1000;       %to convert from mm/s to m/s

open('v_vs_time_fit.mat');
v_vs_time_fit=ans.v_vs_time_fit;
v=v_vs_time_fit.yfit;
%v=v/1000;       %to convert from mm/s to m/s

% generate vel array (=sqrt(v^2)+(v^2))
vel=(sqrt((u.^2)+(v.^2)))

% use central difference method to derive acceleration data from velocities
t=2;            % start at row 1 of WL array

for j=1:(J-2)
    u_prime(j)=(u(t+1)-u(t-1))/(time(t+1)-time(t-1));
    v_prime(j)=(v(t+1)-v(t-1))/(time(t+1)-time(t-1));
    vel_prime(j)=(vel(t+1)-vel(t-1))/(time(t+1)-time(t-1));
t=t+1;
end

u_max=max(u(1:16));
v_max=min(v);
vel_max=max(vel(1:16));
u_v_vel_max=[u_max v_max vel_max]

figure(1)
plot(time,u,'b');
hold on;
plot(time,v,'r');
plot(time,vel,'g');
hold off;
title('Landslide CoM Velocity Time History');
xlabel('time (s)');
ylabel('velocity (mm/s)');
axis([0 7 -500 2000]);
grid on;
legend('u','v','vel',1)

figure(2)
plot(time(2:(J-1)),u_prime,'b');
hold on;
plot(time(2:(J-1)),v_prime,'r');
plot(time(2:(J-1)),vel_prime,'g');
hold off;
title('Landslide CoM Acceleration Time History');
xlabel('time (s)');
ylabel('acceleration (mm/s^2)');
axis([0 7 -2500 2000]);
grid on;
legend('u-prime','v-prime','vel-prime',1)