% Generates bottom profile (boundary condition) time history
% Langford Sue
% 19 April 2006

clear

% Loads landslide Trailing Edge (TE) x-position and time arrays
open('Tx.mat');
Tx=ans.Tx;      % x-position time history of the leading edge of landslide
time=ans.time;  % time corresponding to Lx
[I,J]=size(Tx);

for i=1:J;  % J is the no. of time steps in the TE position time history 
    
    if Tx(i)<=1.297;    % at each time step, calculates Trailing Edge y-pos
        Ty(i)=-(Tx(i)*tan(15*pi/180));  % TE on planar slope
    elseif Tx(i)>=1.807;
        Ty(i)=-0.435;                   % TE on horizontal flume floor
    else                                % TE on transition curve
        Ty(i)=(0.19*((Tx(i)-1.297)^3))+(0.1024*((Tx(i)-1.297)^2))-(.2728*(Tx(i)-1.297))-0.3475;
    end
    
    if Tx(i)<=(1.297-0.5*cos(15*pi/180));   % at each time step of TE pos, calculates corresponding leading edge x- & y-pos
        Lx(i)=Tx(i)+0.500*cos(15*pi/180);   
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
        
    % based on the LE and TE positions, calculates corresponding block CoM    
    Cx(i)=(Tx(i)+((Lx(i)-Tx(i))/2))+(0.013*sin(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
    Cy(i)=(Ly(i)+((Ty(i)-Ly(i))/2))+(0.013*cos(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i)))));
    % based on the LE and TE positions, calculates corresponding origin of ellipse
    Ox(i)=(Tx(i)+((Lx(i)-Tx(i))/2));
    Oy(i)=(Ly(i)+((Ty(i)-Ly(i))/2));
    % based on the LE and TE positions, calculates corresponding angle (degrees) of
    % block to the horizontal
    theta(i)=(atan((Ty(i)-Ly(i))/(Lx(i)-Tx(i))))/(pi/180);
    % check distance between LE and TE
    b(i)=sqrt(((Lx(i)-Tx(i))^2)+((Ly(i)-Ty(i))^2));
    
    % At each time step, create bottom profile at 1mm resolution.    
    Lx(i)=1000.*Lx(i);    % convert to mm, so can round up/down to nearest mm
    Ly(i)=1000.*Ly(i);
    Ox(i)=1000.*Ox(i);
    Oy(i)=1000.*Oy(i);
    Tx(i)=1000.*Tx(i);
    Ty(i)=1000.*Ty(i);
        
    if Tx(i)<=1297-(500*cos(15*pi/180))  % if block is still entirely on planar slope
        for j=1:(301+ceil(Tx(i)))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
            e=((ceil(Tx(i)+1)):1:fix(Lx(i)));   % generated translated semi-ellipse shape
            [K,L]=size(e);
            for k=1:L;
                E_search_min(k)=(Ty(i)-((e(k)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
                E_search_max(k)=(E_search_min(k)+(26/cos(theta(i)*pi/-180))+0.000);
                E_search=[E_search_max(k)];
                string=['(((((',num2str(e(k)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(250^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(k)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(26^2))-1'];
                E(i,k)=fzero(vectorize (inline (char (string))),E_search);
            end      
        j=((301+ceil(Tx(i)))+1);
        for k=1:L   % for transformed half-ellipse
            bottom(i,j)=E(i,k);
            j=j+1;
        end
        for j=(301+ceil(Lx(i))):(301+(1297-1))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
        for j=(301+1297):(301+(1807-1))   % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
        for j=(301+1807):(301+10200)   % horizontal flume floor
            bottom(i,j)=-435;
        end
    elseif Tx(i)>1297-(500*cos(15*pi/180)) & Tx(i)<=1297;   % if block is part on the planar slope, part on the transition curve
        for j=1:(301+ceil(Tx(i)))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
            e=((ceil(Tx(i)+1)):1:fix(Lx(i)));   % generated translated semi-ellipse shape
            [K,L]=size(e);
            for k=1:L;
                E_search_min(k)=(Ty(i)-((e(k)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
                E_search_max(k)=(E_search_min(k)+(26/cos(theta(i)*pi/-180))+0.000);
                E_search=[E_search_max(k)];
                string=['(((((',num2str(e(k)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(250^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(k)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(26^2))-1'];
                E(i,k)=fzero(vectorize (inline (char (string))),E_search);
            end      
        j=((301+ceil(Tx(i)))+1);
        for k=1:L   % for transformed half-ellipse
            bottom(i,j)=E(i,k);
            j=j+1;
        end
        for j=(301+ceil(Lx(i))):(301+(1807-1))   % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
        for j=(301+1807):(301+10200)   % horizontal flume floor
            bottom(i,j)=-435;
        end
    elseif Tx(i)>1297 & Tx(i)<=1314;   % if block is entirely on the transition curve
        for j=1:(301+(1297-1))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
        for j=(301+1297):(301+ceil(Tx(i)))     % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
            e=((ceil(Tx(i)+1)):1:fix(Lx(i)));   % generated translated semi-ellipse shape
            [K,L]=size(e);
            for k=1:L;
                E_search_min(k)=(Ty(i)-((e(k)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
                E_search_max(k)=(E_search_min(k)+(26/cos(theta(i)*pi/-180))+0.000);
                E_search=[E_search_max(k)];
                string=['(((((',num2str(e(k)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(250^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(k)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(26^2))-1'];
                E(i,k)=fzero(vectorize (inline (char (string))),E_search);
            end      
        j=((301+ceil(Tx(i)))+1);
        for k=1:L   % for transformed half-ellipse
            bottom(i,j)=E(i,k);
            j=j+1;
        end
        for j=(301+ceil(Lx(i))):(301+(1807-1))   % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
        for j=(301+1807):(301+10200)   % horizontal flume floor
            bottom(i,j)=-435;
        end
    elseif Tx(i)>1314 & Tx(i)<=1807;   % if block is partly on the transition curve, partly on the flume floor
        for j=1:(301+(1297-1))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
        for j=(301+1297):(301+ceil(Tx(i)))     % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
            e=((ceil(Tx(i)+1)):1:fix(Lx(i)));   % generated translated semi-ellipse shape
            [K,L]=size(e);
            for k=1:L;
                E_search_min(k)=(Ty(i)-((e(k)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
                E_search_max(k)=(E_search_min(k)+(26/cos(theta(i)*pi/-180))+0.000);
                E_search=[E_search_max(k)];
                string=['(((((',num2str(e(k)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(250^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(k)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(26^2))-1'];
                E(i,k)=fzero(vectorize (inline (char (string))),E_search);
            end          
        j=((301+ceil(Tx(i)))+1);
        for k=1:L   % for transformed half-ellipse
            bottom(i,j)=E(i,k);
            j=j+1;
        end
        for j=(301+ceil(Lx(i))):(301+10200)   % horizontal flume floor
            bottom(i,j)=-435;
        end
    elseif Tx(i)>1807 & Tx(i)<=10200;   % if block is entirely on the horizontal flume floor
        for j=1:(301+(1297-1))   % planar slope
            bottom(i,j)=-(j-301)*tan(15*pi/180);
        end
        for j=(301+1297):(301+(1807-1))   % transition curve
            bottom(i,j)=1000*((0.19*((((j-301)/1000)-1.297)^3))+(0.1024*((((j-301)/1000)-1.297)^2))-(.2728*(((j-301)/1000)-1.297))-0.3475);
        end
        for j=(301+1807):(301+ceil(Tx(i)))   % horizontal flume floor
            bottom(i,j)=-435;
        end
            e=((ceil(Tx(i)+1)):1:fix(Lx(i)));   % generated translated semi-ellipse shape
            [K,L]=size(e);
            for k=1:L;
                E_search_min(k)=(Ty(i)-((e(k)-Tx(i))*(Ty(i)-Ly(i))/(Lx(i)-Tx(i))));
                E_search_max(k)=(E_search_min(k)+(26/cos(theta(i)*pi/-180))+0.000);
                E_search=[E_search_max(k)];    
                string=['(((((',num2str(e(k)),'-',num2str(Ox(i)),')*cos(',num2str(theta(i)*pi/-180),'))+((y-',num2str(Oy(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(250^2))+(((((y-',num2str(Oy(i)),')*cos(',num2str(theta(i)*pi/-180),'))-((',num2str(e(k)),'-',num2str(Ox(i)),')*sin(',num2str(theta(i)*pi/-180),')))^2)/(26^2))-1'];
                E(i,k)=fzero(vectorize (inline (char (string))),E_search);
            end              
        j=((301+ceil(Tx(i)))+1);
        for k=1:L   % for transformed half-ellipse
            bottom(i,j)=E(i,k);
            j=j+1;
        end
        for j=(301+ceil(Lx(i))):(301+10200)   % horizontal flume floor
            bottom(i,j)=-435;
        end
    end
    i=i+1
end 

x=-300:1:10200; %creates x-position array

figure(1)
set(axes('FontSize', 12));
colormap(jet)
surf(x,time,bottom), shading interp;     % plots bottom profile vs space and time
title('Bottom Profile Time History');
xlabel('distance from shore (mm)');
ylabel('time (s)');
zlabel('d (mm)');
axis([-300 10200 0 6 -500 100]);    % Fix up axes to fit image