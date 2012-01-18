%% Reader cylindrical coordinate 2D

clear all
A0 = importdata('out0.txt');
A1 = importdata('outHankelError.txt');

%%

time=A1(:,1);
ErrorBoth=A1(:,2);
ErrorHankel=A1(:,3);
ErrorFFT=A1(:,4);

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');


plot(time,ErrorBoth,'LineWidth',3)
hold on
plot(time,ErrorHankel,'ro','LineWidth',1)
hold on
plot(time,ErrorFFT,'mx','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Error','fontsize',12)
title('Error in normalization','fontsize',16)
grid on
hleg1 = legend('Both','Hankel','FFT');

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');

plot(time,log10(abs(ErrorBoth)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorHankel)),'ro','LineWidth',1)
hold on
plot(time,log10(abs(ErrorFFT)),'mx','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Logaritmic error in normalization','fontsize',16)
grid on
hleg1 = legend('Both','Hankel','FFT');