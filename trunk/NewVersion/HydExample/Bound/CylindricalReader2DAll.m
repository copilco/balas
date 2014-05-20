%% Reader cylindrical coordinate 2D

clear all

ax  = importdata('axes.txt');
%qax = importdata('qaxes.txt');
A0 = importdata('out0.txt');
A2 = importdata('out2.txt');




%% Grid parameters


skiper=1;


nr = 300/skiper; %A3(1);
nz = 600/skiper; %A3(2);


%kr = qax(1:skiper:nr);
%kz = qax(nr+1:skiper:nz+nr);

Nr = nr;%length(kr) ;
Nz = nz;%length(kz);



snap  = length(A0)/Nr/Nz;%1;%
Nsnap = snap-1;%floor(Ntime/snap);




r = ax(1:skiper:skiper*nr);
z = ax(skiper*nr+1:skiper:skiper*(nz+nr));




[Z R]   = meshgrid(z,r);
%vpot    = reshape(A2,Nz,Nr);



%[KZ KR] = meshgrid(kz,kr);



%% Movie time interaction
% aviobj = avifile('EWP.avi');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.8],...
    'Color','w');

ksnap = 20;

xmin = -40;
xmax =  40;
ymin =  0;
ymax =  40;


krmax = 4;
krmin = 0;
kzmax = 4;
kzmin = -4;
wx=0.057;
Ip=0.5;
p0=sqrt(2*(wx-Ip));


for j=1:Nsnap
   
    
    SQ_PHI       = reshape(A0(1+Nr*Nz*(j-1):Nr*Nz*j),Nz,Nr);
    
    %subplot(1,2,1)
    surfc(Z,R,log10(SQ_PHI'+1e-9),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    xlabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    ylabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    
    
    
    
    title(['ktime= ',num2str(j),' Full Distribution ' ],'fontsize',16)
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
        'fontsize',16)%,'fontweight','b');    
    
    %caxis([-12 0])
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)
    

    
%     subplot(1,2,2)
%     surfc(Z,R,vpot.',...
%         'FaceColor','interp',...
%         'EdgeColor','none')
%     
%     view(2)
%     axis tight
%     
%     h=gca;    
%     set(h,'fontsize',16)
%     
%     xlabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
%     ylabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
%     
%     
%     
%     
%     %title(['ktime= ',num2str(ksnap*j),' Full Distribution ' ],'fontsize',16)
%     h = colorbar('location','EastOutside');
%     set(get(h,'YLabel'),'String',' V(z,\rho)  ',...
%         'fontsize',16)%,'fontweight','b');    
%     
%     %caxis([-12 0])
%     
%     xlim([xmin xmax])
%     ylim([ymin ymax])
%     
%     h=gca;    
%     set(h,'fontsize',16)    
    
    
    

    namefile=['SnapshotPos',num2str(j),'.jpg'];
    %saveas(gcf,namefile, 'jpg')
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [16 10]);
    set(gcf, 'PaperPosition', [0 0 16 10]);
    %name=['soft_eff' num2str(soft_eff(num),2),'.jpg'];
    print('-djpeg','-cmyk',namefile);     
    
    pause(0.01)
end




%%
figure
plot(z,sum(SQ_PHI(:,1:40),2))


%%
figure
plot(z,vpot(:,1))

%%
break



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


%% 
static_dipole = sum(sum((SQ_PHI.').*Z))


%%