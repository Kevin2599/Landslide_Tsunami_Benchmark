ChartTitle='Full Wave Field - SG3-IS5';            %3D plot title

% Loads Position, Time, and Water Level variables from saved combined .mat
% file
open('PosAll_1_31.mat');
Zall_Xmod=ans.Zall_Xmod./1000;
Zall_WL=ans.Zall_WL./1000;
Zall_T=ans.Zall_T;

% Create new meshgrid with smaller time intervals and larger spatial
% intervals and interpolate WL on this new grid to remove 'hairiness'.
% Cubic used to interpolate
[Position,Time] = meshgrid(0:0.00399:10.100,0:.006667:8.0004);
WaterLevel = interp2(Zall_Xmod,Zall_T,Zall_WL,Position,Time,'cubic');

% Reads in x position time history.  Slider block CoM location
open('x_vs_time_fit.mat');
x_vs_time_fit=ans.x_vs_time_fit;
time=x_vs_time_fit.xi;
x=x_vs_time_fit.yfit/1000;

% creates false 3rd dimesion to CoM position data to move it above the WL
% plot
for i=1:120;
    y_false(i)=40;
end
% pads CoM x-position with last position value
for i=91:120;
    x(i)=x(90);
end
% pads time array with full 8sec
for i=1:120;
    time(i)=(i-1)*(1/15);
end

% Plots 2D WL profile
figure(1)
set(axes('FontSize', 14));
colormap(gray)
%contourf(Position,Time,WaterLevel,14);     % plots combined zone water level profile time history
mesh(Position,Time,WaterLevel);
colorbar;
shading interp;

hold on;
plot3(x,time,y_false,'k o-');     % plots block CoM position on graph
%plot(maxpos,Zall_T,'b');
%plot(minpos,Zall_T,'r');
hold off;

hold on;
run 'CharacteristicCurve_GroupVel';
hold off;

Title=title(ChartTitle);
xlabel('x (m)');
ylabel('time (s)');
%zlabel('Water level (mm)');
%legend('Block CoM position',2)
legend('landslide CoM position','group velocity (wavelength=0.15m)','phase velocity (wavelength=4m)',4)
%legend('landslide CoM position','phase characteristic lambda=0.2m','group characteristic lambda=0.2m',...
%    'phase characteristic lambda=0.5m','group characteristic lambda=0.5m',...
%    'phase characteristic lambda=1m','group characteristic lambda=1m',...
%    'phase characteristic lambda=2m','group characteristic lambda=2m',...
%    'phase characteristic lambda=4m','group characteristic lambda=4m',...
%    'phase characteristic lambda=6m','group characteristic lambda=6m',4)
%legend('landslide CoM position','phase characteristic lambda=0.15m','group characteristic lambda=0.15m',...
%    'phase characteristic lambda=4m','group characteristic lambda=4m',...
%    'phase characteristic lambda=6m','group characteristic lambda=6m',4)
axis([-0.100 10.200 0 8]);    % Fix up axes to fit image

clear i m n