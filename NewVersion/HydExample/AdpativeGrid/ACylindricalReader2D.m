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
Aa = importdata('outa.txt'); 
A7 = importdata('out7.txt');
A9 = importdata('out9.txt');
A10 = importdata('out10.txt');


jparam = 1;

%% Grid parameters
nr      = 2500/4;%1700; %A3(1);
nz      = 4500/4;%3000; %A3(2);

snaper1 =1;
kr      = qax(1:snaper1:nr);
kz      = qax(nr+1:snaper1:nz+nr);

Nr      = length(kr) ;
Nz      = length(kz);

snap    = length(A3)/Nr/Nz;%1;%
Nsnap   = snap-1;%floor(Ntime/snap);

r       = ax(1:1:nr);
z       = ax(nr+1:1:nz+nr);

[Z R]   = meshgrid(z,r);


[KZ KR] = meshgrid(kz,kr);

dt      = 0.05;

t       = Aa(:,1);
Exuv    = Aa(:,4);
ksnap   = floor(length(t)/snap);


%%%
%break

%% Movie time interaction
% aviobj = avifile('EWP.avi');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.85 scrsz(4)*0.8],...
    'Color','w');
krmax = 2;
krmin = 0.1;
kzmax = 2;
kzmin = -2;
wx=0.057;
Ip=0.5;
p0=sqrt(2*(wx-Ip));

theta=0: (2*pi/(nz-1)): 2*pi;


xmin = min(z);
xmax =  max(z);
ymin =  min(r);
ymax =  max(r);
zerowave = 1e-15;
qzero    = zerowave;
for j=1:Nsnap
   
    
    %clf
    
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/1.25],...
%         'Color','w');
    
    PHI  = reshape(A3(1+Nr*Nz*(j-1):Nr*Nz*j),Nz,Nr);
    PHI_Mask  = reshape(A4(1+Nr*Nz*(j-1):Nr*Nz*j),Nz,Nr);
    
    
    subplot(3,4,[1 2 3 4] )
    surfc(Z,R,log10(PHI'+zerowave),...
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
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)
    

    
    
    subplot(3,4,[5 6 7 8])
    surfc(Z,R,log10(PHI_Mask'+zerowave),...
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
    
    caxis(log10([zerowave 1e-1]))
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)    
    
    
    subplot(3,4, 9)
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





%scrsz = get(0,'ScreenSize');
%figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.8],...
%    'Color','w');



    
    %clf
    
    %     scrsz = get(0,'ScreenSize');
    %     figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/1.25],...
    %         'Color','w');
    
    subplot(3,4,[10 11 12])
    PHIQ = reshape(A5(1+Nr*Nz*(j-1):Nr*Nz*j),Nz,Nr);
    
    
    surf(KZ,KR,log10(abs(PHIQ')+qzero),'FaceColor','interp','EdgeColor','none');
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
    



    ylabel('k\rho (a.u.) ','fontsize',16) %,'fontweight','b')
    xlabel('kz (a.u.) ','fontsize',16)    %,'fontweight','b')
    %caxis(log10([qzero 1e-5]))
    
    title('ElectronicWavePacket in momentum space','fontsize',16)
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
        'fontsize',16)%,'fontweight','b');
    
    %caxis([-10 -6.5])
    
    h=gca;
    set(h,'fontsize',16)    
    
    nameofjpg = sprintf('ATITest0_%d.jpg',j);
    %saveas(gcf,nameofjpg,'jpg')
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','inches');    
    set(gcf,'PaperSize',[16 10]);        
    set(gcf,'PaperPosition',[0 0 16 10]);    
    print('-djpeg','-cmyk',nameofjpg);
    
    
    pause(0.01)
end


!mkdir Figure
!mv *.jpg


%% Final momentum distribution %% 

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.2 scrsz(4)*0.8],...
    'Color','w');

plot(kz,PHIQ(:,1),'LineWidth',1.5)
xlim([kzmin kzmax])

grid on

xlabel('Kz (a.u.)','fontsize',16)
ylabel('Momentum distribution','fontsize',16)
title('Momentum distribution in z axis','fontsize',16)



%% Dipole Acceleration and Spectrum depicted %% 


I0      = A9(1)/3.5e16;
w0      = 0.057;

Up      = I0/4/w0^2;
Ip      = abs(A9(7));

nmax    = floor((Ip+3.17*Up)/w0);

tp      = A7(:,1);
dtp     = tp(2)-tp(1);
Ntp     = length(tp);

wmaxp   = pi/dtp;
dwp     = 2*wmaxp/Ntp;

wp      = (-wmaxp:dwp:(wmaxp-dwp)).';

adip_time = A7(:,3);
adip_wp   = dtp*fftshift(fft(ifftshift(adip_time)));


szero           = 1e-12;
adip_spectrum   = abs(adip_wp).^2+szero;



%E0 = ;

figure(1)
plot(A7(:,1),A7(:,2),'r','linewidth',2)
hold on
plot(A7(:,1),A7(:,3),'g','linewidth',2)
xlabel('time (a.u.) ','fontsize',16)
set(gca,'fontsize',16)
ylabel('a_z (a.u.) ','fontsize',16)

axis tight


nameofjpg = sprintf('AccelDipole01_%d.jpg',jparam);
%saveas(gcf,nameofjpg,'jpg')
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperSize',[9 4]);
%set(gcf,'PaperFont',16);
set(gcf,'PaperPosition',[0 0 9 4]);    
print('-djpeg','-cmyk',nameofjpg);
%print('-djpeg', '-f1', '-r300', nameofjpg)


%% Spectrum

figure
plot(wp/w0, log10(adip_spectrum),'linewidth',2)
hold on
%line([1 1]*nmax,log10([szero 1e-4]),'Color','r','linestyle','--')
line([1 1]*nmax,[-12 3],'Color','r','linestyle','--')
xlabel('Harmonic-Order ','fontsize',15)
ylabel('log_{10}(|a_z(\omega)|^2) ','fontsize',16)

legend(['I_0= ',num2str(I0*3.5e16,'%1.2e'), 'W/cm^2 '],4)

xlim([0 25])
ylim([-7 3])
set(gca,'fontsize',16)


nameofjpg = sprintf('HHGSpectraDipole0_%d.jpg',jparam);
%saveas(gcf,nameofjpg,'jpg')
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');    
set(gcf,'PaperSize',[16 10]);        
set(gcf,'PaperPosition',[0 0 16 10]);    
print('-djpeg','-cmyk',nameofjpg);



DATA=[tp A7(:,2) A7(:,3) wp real(adip_wp) imag(adip_wp) abs(adip_wp).^2];


filename = ['Balas1s2p_0' num2str(jparam) '.txt'];

save(filename,'DATA','-ASCII');

!cp Balas1s2p* ../ComparionIntensiy

%% Population Vs time %% 

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*.5 scrsz(4)*0.48],...
    'Color','w');

plot(A10(:,1),A10(:,2)/max(abs(A10(:,2)))*0.35+.5,'m','LineWidth',1.5)
%xlim([kzmin kzmax])
hold on
plot(A10(:,1),A10(:,3),'r','LineWidth',1.5)
plot(A10(:,1),A10(:,4),'b','LineWidth',1.5)

grid on
set(gca,'fontsize',16)
xlabel('time (a.u.)','fontsize',16)
ylabel('Population ','fontsize',16)
title('Populations ','fontsize',16)

%%