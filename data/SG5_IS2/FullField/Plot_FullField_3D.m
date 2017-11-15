ChartTitle='Full Wave Field - SG5-IS2';            %3D plot title

% Loads Position, Time, and Water Level variables from saved combined .mat
% file
open('PosAll_1_31.mat');
Zall_Xmod=ans.Zall_Xmod;
Zall_WL=ans.Zall_WL;
Zall_T=ans.Zall_T;

% Create new meshgrid with smaller time intervals and larger spatial
% intervals and interpolate WL on this new grid to remove 'hairiness'.
% Cubic used to interpolate
[Position,Time] = meshgrid(0:3.99:10100,0:.006667:8.0004);
WaterLevel = interp2(Zall_Xmod,Zall_T,Zall_WL,Position,Time,'cubic');

% Reads in x position time history.  Slider block CoM location
open('x_vs_time_fit.mat');
x_vs_time_fit=ans.x_vs_time_fit;
time=x_vs_time_fit.xi;
x=x_vs_time_fit.yfit;

% Plots 3D WL profile
figure(1)
set(axes('FontSize', 12));
colormap(jet)
surf(Position,Time,WaterLevel), shading interp;     % plots combined zone water level profile time history

hold on;
plot(x,time,'k o-');     % plots block CoM position on graph
hold off;

Title=title(ChartTitle);
xlabel('Distance from shore (mm)');
ylabel('Time (s)');
zlabel('Water level (mm)');
legend('Block CoM position',2)
axis([-100 10200 0 8 -30 30]);    % Fix up axes to fit image
