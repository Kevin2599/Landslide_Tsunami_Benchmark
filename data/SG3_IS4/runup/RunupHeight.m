ChartTitle='Wave Run-up Time History';            %3D plot title

% Loads Position, Time, and Water Level variables from saved combined .mat
% file
open('SG3_IS4_r1_12pt582_200.mat');
WL=ans.WL;
T=ans.T;
maxima=ans.maxima;

WL=WL.*sin(15*pi/180);

% Plots wave run-up time history
figure(1)
set(axes('FontSize', 12));
plot(T,WL./1000,'ko:'); 
Title=title(ChartTitle);
xlabel('time (s)');
ylabel('wave run-up height (m)');
%legend('100','150','200','240',2)
grid on;
axis([0 8 -0.015 0.010]);    % Fix up axes to fit image