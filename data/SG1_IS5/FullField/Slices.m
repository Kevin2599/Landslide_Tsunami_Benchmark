
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

