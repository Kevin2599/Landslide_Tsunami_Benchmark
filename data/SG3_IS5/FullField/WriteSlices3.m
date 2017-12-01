
%--------------------------------------------------------------------------
load('PosAll_1_31.mat');
%Zall_Xmod=ans.Zall_Xmod;
%Zall_WL=ans.Zall_WL/1000;
%Zall_T=ans.Zall_T;

Zall_Xmod=Zall_Xmod./1000;
Zall_WL=Zall_WL/1000;

y1=Zall_WL(25,(120:end));
y2=Zall_WL(40,(120:end));
y3=Zall_WL(55,(120:end));
y4=Zall_WL(70,(120:end));
y5=Zall_WL(85,(120:end));
y6=Zall_WL(100,(120:end));

x=Zall_Xmod(120:end);

% create a matrix y, with 7 rows
y = [x(:), y1(:), y2(:), y3(:), y4(:), y5(:), y6(:)];

save('profile.dat','y','-ascii')


y1=Zall_WL((16:end),1373);
y2=Zall_WL((16:end),3879);
y3=Zall_WL((16:end),6386);
y4=Zall_WL((16:end),8892);
y5=Zall_WL((16:end),11398);
y6=Zall_WL((16:end),13904);

t=Zall_T(16:end);

% create a matrix y, with 7 rows
yy = [t(:), y1(:), y2(:), y3(:), y4(:), y5(:), y6(:)];
save('timeseries.dat','yy','-ascii')
 

%--------------------------------------------------------------------------

