%% Reader cylindrical coordinate 2D

clear all
A0 = importdata('out0.txt');
A1 = importdata('outCrankUniError.txt');

%%

time=A1(:,1);
ErrorBoth=A1(:,2);
ErrorRHO=A1(:,3);
ErrorZ=A1(:,4);

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.9],'Color','w');

plot(time,ErrorBoth,'LineWidth',3)
hold on
plot(time,ErrorRHO,'ro','LineWidth',1)
hold on
plot(time,ErrorZ,'mx','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Error','fontsize',12)
title('Error in normalization','fontsize',16)
grid on
hleg1 = legend('Both','rho','z');


%%
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.9],'Color','w');

plot(time,log10(abs(ErrorBoth)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorRHO)),'ro','LineWidth',1)
hold on
plot(time,log10(abs(ErrorZ)),'mx','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Logaritmic error in normalization','fontsize',16)
grid on
hleg1 = legend('Both','rho','z');