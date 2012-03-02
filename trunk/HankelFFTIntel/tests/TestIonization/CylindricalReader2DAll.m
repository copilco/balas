%% Reader cylindrical coordinate 2D

clear all

ax = importdata('axes.txt');
qax = importdata('qaxes.txt');
%A0 = importdata('out0.txt');
%A1 = importdata('out1.txt');
%A2 = importdata('out2.txt');
A3 = importdata('out3.txt');
A4 = importdata('out4.txt'); 
A5 = importdata('out5.txt'); 
A6 = importdata('out6.txt'); 

%% Grid parameters
nr = 100%500; %A3(1);
nz = 200%1000; %A3(2);


snap  = length(A3)/nr/nz;%1;%
Nsnap = snap-1;%floor(Ntime/snap);

r = ax(1:nr);
z = ax(nr+1:nz+nr);

[Z R] = meshgrid(z,r);

kr = qax(1:nr);
kz = qax(nr+1:nz+nr);

[KZ KR] = meshgrid(kz,kr);

dt = 0.05;

t     = A6(:,1);
Exuv  = A6(:,2);
ksnap = floor(length(t)/snap);



%% Movie time interaction
% aviobj = avifile('EWP.avi');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.8],...
    'Color','w');


xmin = -50;
xmax =  50;
ymin =  0;
ymax =  60;


for j=1:4:Nsnap+2
    
    %clf
    
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/1.25],...
%         'Color','w');
    
    PHI  = reshape(A3(1+nr*nz*(j-1):nr*nz*j),nz,nr);
    PHI_Mask  = reshape(A4(1+nr*nz*(j-1):nr*nz*j),nz,nr);
    
    
    subplot(3,4,[1 2 3 4] )
    surfc(Z,R,log10(PHI'+1e-12),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    %xlabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    ylabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    
    
    
    
    title(['ktime= ',num2str(ksnap*j),' Full Distribution ' ],'fontsize',16)
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
        'fontsize',16)%,'fontweight','b');    
    
    %caxis([-12 0])
    
    %xlim([xmin xmax])
    %ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)
    

    
    
    subplot(3,4,[5 6 7 8])
    surfc(Z,R,log10(PHI_Mask'+1e-12),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    xlabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    ylabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    
    
    
    
    title('Masking ','fontsize',16)
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
        'fontsize',16)%,'fontweight','b');    
    
    caxis([-12 0])
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)    
    
    
    subplot(3,4,[ 9])
    plot(t,Exuv/max(Exuv),'--b','linewidth',1)
    hold on
    plot(t(1:ksnap*j) ,Exuv(1:ksnap*j)/max(Exuv),'Color',[0.5 0 1],'linewidth',6)
    
    
    xlabel('t (a.u.) ','fontsize',16)%   
    ylabel('Exuv (arb. u.) ','fontsize',16)%
    
    
    h=gca;    
    set(h,'fontsize',16)    
    axis tight    
    %F = getframe(fig);
    %aviobj = addframe(aviobj,F); 
    
    %pause(0.01)



krmax = 2;
krmin = 0;
kzmax = 2;
kzmin = -2;
wx=0.057;
Ip=0.5;
p0=sqrt(2*(wx-Ip));

theta=0: (2*pi/(nz-1)): 2*pi;


scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.8],...
%    'Color','w');



    
    %clf
    
    %     scrsz = get(0,'ScreenSize');
    %     figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/1.25],...
    %         'Color','w');
    
    subplot(3,4,[10 11 12])
    PHIQ = reshape(A5(1+nr*nz*(j-1):nr*nz*j),nz,nr);
    
    
    surf(KZ,KR,log10(abs(PHIQ')+1e-11),'FaceColor','interp','EdgeColor','none');
    view(2)
    
    %view([90 -90])
    axis tight
    %daspect([1,1,1])
    
    
    
   %  hold on
%      for h=9:20
%          p0=sqrt(2*(h*wx-Ip));
%           plot3(p0*cos(theta), p0*sin(theta),ones(1,nz)*(-9),'linewidth',3,'Color','w' )
%      end
          ylim([krmin krmax])
     xlim([kzmin kzmax])
    

    h=gca;
    set(h,'fontsize',16)

    xlabel('k\rho (a.u.) ','fontsize',16)%,'fontweight','b')
    ylabel('kz (a.u.) ','fontsize',16)%,'fontweight','b')

    title('ElectronicWavePacket in momentum space','fontsize',16)
    h = colorbar('location','EastOutside');
    %set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
    %    'fontsize',16)%,'fontweight','b');
    %caxis([-10 -6.5])

     nameofjpg = sprintf('ATITest0_%d.jpg',j);
    saveas(gcf,nameofjpg,'jpg')
    
    pause(0.2)
end


%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.2 scrsz(4)*0.8],...
    'Color','w');

plot(kz,PHIQ(:,1),'LineWidth',1.5)
xlim([kzmin kzmax])

grid on

xlabel('Kz (a.u.)','fontsize',16)
ylabel('Momentum distribution','fontsize',16)
title('Momentum distribution in z axis','fontsize',16)
