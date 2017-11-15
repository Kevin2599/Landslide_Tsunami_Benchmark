% Generates bottom profile (boundary condition) time history

clear

% Loads x and y position, and u and v velocity, and time
open('x_vs_time_fit_model.mat');
x_vs_time_fit=ans.x_vs_time_fit_model;
time=x_vs_time_fit.xi;
x=x_vs_time_fit.yfit;
x=x/1000;       %to convert from mm to m
[I,J]=size(x);

open('y_vs_time_fit_model.mat');
y_vs_time_fit=ans.y_vs_time_fit_model;
y=y_vs_time_fit.yfit;
y=y/1000;       %to convert from mm to m


for i=1:I;  % at each block position in each frame, 
    K=(0.2800:-.0005:0.2100);   %trials trailing edge of block from 0.28m-0.21m behind measured block CoM 
    for k=1:141;
        Tx(i)=x(i)-K(k);
            if Tx(i)<=1.297;    % at each trial for trailing edge x-pos, calculates TE y-pos
                Ty(i)=-(Tx(i)*tan(15*pi/180));
            elseif Tx(i)>=1.807;
                Ty(i)=-0.435;
            else 
                Ty(i)=(0.19*((Tx(i)-1.297)^3))+(0.1024*((Tx(i)-1.297)^2))-(.2728*(Tx(i)-1.297))-0.3475;
            end

                if Tx(i)<=(1.297-0.5*cos(15*pi/180));
                    Lx(i)=Tx(i)+0.500*cos(15*pi/180);   % for each trial for TE pos, calculates corresponding leading edge x- & y-pos
                    Ly(i)=Ty(i)-0.500*sin(15*pi/180);
                elseif Tx(i)>(1.297-0.5*cos(15*pi/180)) & Tx(i)<1.314;
                    string=['((x-',num2str(Tx(i)),')^2)+((',num2str(Ty(i)),'-((0.19*((x-1.297)^3))+(0.1024*((x-1.297)^2))-(0.2728*(x-1.297))-0.3475))^2)-(0.500^2)'];
                    Lx(i)=fzero(vectorize (inline (char (string))),1.8);
                    Ly(i)=(0.19*((Lx(i)-1.297)^3))+(0.1024*((Lx(i)-1.297)^2))-(0.2728*(Lx(i)-1.297))-0.3475;
                elseif Tx(i)>=1.314 & Tx(i)<1.807 
                    Ly(i)=-0.435;
                    Lx(i)=(sqrt((0.500^2)-((Ty(i)-Ly(i))^2)))+Tx(i);
                elseif Tx(i)>=1.807;
                    Lx(i)=Tx(i)+0.500;
                    Ly(i)=-0.435;
                else
                    Lx(i)=0;
                    Ly(i)=0;
                end
            
            % at each trial of TE position, calculates corresponding block centre    
            Cx(i)=(Tx(i)+((Lx(i)-Tx(i))/2))+(0.013*sin(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
            Cy(i)=(Ly(i)+((Ty(i)-Ly(i))/2))+(0.013*cos(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
         
        % at each trial of TE position, compares calculated CoM position with measured CoM position    
        %delta(k)=sqrt(((Cx(i)-x(i))^2)+((Cy(i)-y(i))^2));   % minimise absolute distance difference
        delta(k)=sqrt(((Cy(i)-y(i))^2));                   % minimise vertical distance difference
        [delta_min(i),Kmin(i)]=min(delta);  % finds at what trial for TE position provides the CoM closest to the measured CoM
        
        % calculate Tx,Ty,Lx,Ly,Cx,Cy at the optimal TE trial location
        Tx(i)=x(i)-K(Kmin(i));
            if Tx(i)<=1.297;    
                Ty(i)=-(Tx(i)*tan(15*pi/180));
            elseif Tx(i)>=1.807;
                Ty(i)=-0.435;
            else 
                Ty(i)=(0.19*((Tx(i)-1.297)^3))+(0.1024*((Tx(i)-1.297)^2))-(.2728*(Tx(i)-1.297))-0.3475;
            end

                if Tx(i)<=(1.297-0.5*cos(15*pi/180));
                    Lx(i)=Tx(i)+0.500*cos(15*pi/180);   % for each trial for TE pos, calculates corresponding leading edge x- & y-pos
                    Ly(i)=Ty(i)-0.500*sin(15*pi/180);
                elseif Tx(i)>(1.297-0.5*cos(15*pi/180)) & Tx(i)<1.314;
                    string=['((x-',num2str(Tx(i)),')^2)+((',num2str(Ty(i)),'-((0.19*((x-1.297)^3))+(0.1024*((x-1.297)^2))-(0.2728*(x-1.297))-0.3475))^2)-(0.500^2)'];
                    Lx(i)=fzero(vectorize (inline (char (string))),1.8);
                    Ly(i)=(0.19*((Lx(i)-1.297)^3))+(0.1024*((Lx(i)-1.297)^2))-(0.2728*(Lx(i)-1.297))-0.3475;
                elseif Tx(i)>=1.314 & Tx(i)<1.807 
                    Ly(i)=-0.435;
                    Lx(i)=(sqrt((0.500^2)-((Ty(i)-Ly(i))^2)))+Tx(i);
                elseif Tx(i)>=1.807;
                    Lx(i)=Tx(i)+0.500;
                    Ly(i)=-0.435;
                else
                    Lx(i)=0;
                    Ly(i)=0;
                end
        
            % at the optimum TE position, calculates corresponding block CoM    
            Cx(i)=(Tx(i)+((Lx(i)-Tx(i))/2))+(0.013*sin(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
            Cy(i)=(Ly(i)+((Ty(i)-Ly(i))/2))+(0.013*cos(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
            % at the optimum TE position, calculates corresponding origin of ellipse
            Ox(i)=(Tx(i)+((Lx(i)-Tx(i))/2));
            Oy(i)=(Ly(i)+((Ty(i)-Ly(i))/2));
            % at the optimum TE position, calculates corresponding angle (degrees) of
            % block to the horizontal
            theta(i)=(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i))))/(pi/180);   
        end 
    
    % at each time step (i), calculate the ellipse profile    
    %e=(Tx(i):0.001:Lx(i));
    %[K,L]=size(e);

    %for j=1:L;
    %    E_search_min(j)=(Ty(i)-((e(j)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
    %    E_search_max(j)=(E_search_min(j)+(0.026/cos(theta(i)*pi/-180))+0.000);
    %    E_search=[E_search_max(j)];
   % 
   %     string=['(((((',num2str(e(j)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(0.25^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(j)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(0.026^2))-1'];
   %     E(i,j)=fzero(vectorize (inline (char (string))),E_search);
   % end
        
    i=i+1
end

clear I J K Kmin ans delta i k string x_vs_time_fit y_vs_time_fit

figure(1)
plot(x,time,'b');
hold on;
plot(Cx,time,'r');
hold off;
title('Comparison of Calculated CoM x-position');
xlabel('distance from shore (m)');
ylabel('time (s)');
axis([0 6 0 6]);
grid on;
legend('measured','calculated',1)

figure(2)
plot(y,time,'b');
hold on;
plot(Cy,time,'r');
hold off;
title('Comparison of Calculated CoM y-position');
xlabel('distance from still water level (m)');
ylabel('time (s)');
axis([-0.5 0 0 6]);
grid on;
legend('measured','calculated',1)

figure(3)
plot(x,y,'b');
hold on;
plot(Cx,Cy,'r');
hold off;
title('Comparison of Calculated CoM Position');
xlabel('distance from shore (m)');
ylabel('distance from still water level (m)');
axis([0 6 -0.5 0]);
grid on;
legend('measured','calculated',1)

figure(4)
plot(time,delta_min,'b');
title('Absolute Distance Between Measured and Calculated CoM Position');
xlabel('time (s)');
ylabel('absolute distance between measured and calculated CoM position (m)');
axis([0 6 0 0.02]);
grid on;

figure(5)
plot(time,theta,'b');
title('Block Angle Relative to Horizontal');
xlabel('time (s)');
ylabel('block angle relative to horizontal (degrees)');
axis([0 6 0 16]);
grid on;