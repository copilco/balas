%% Reader cylindrical coordinate 2D

clear all
%A0 = importdata('out0.txt');
%A1 = importdata('outH2U.txt');
A2 = importdata('outErrorH2U.txt');
%A3 = importdata('outU2H.txt');
A4 = importdata('outErrorU2H.txt');

%%

time=A2(:,1);
%H2U=A1(:,2)-A1(:,3);
%U2H=A3(:,2)-A3(:,3);

ErrorH2Uh=A2(:,2);
ErrorH2Uu=A2(:,3);

ErrorU2Hh=A4(:,2);
ErrorU2Hu=A4(:,3);

%% Plot the error in normalization H2U

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');

plot(time,log10(abs(ErrorH2Uh)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorH2Uu)),'ro','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Comparision in error normalization. H2U','fontsize',16)
grid on
hleg1 = legend('Hankel','Hankel Interp');


%% Plot the error in normalization H2U

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');

plot(time,log10(abs(ErrorU2Hh)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorU2Hu)),'ro','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Comparision in error normalization. U2H','fontsize',16)
grid on
hleg1 = legend('Uniform','Uniform Interp');
