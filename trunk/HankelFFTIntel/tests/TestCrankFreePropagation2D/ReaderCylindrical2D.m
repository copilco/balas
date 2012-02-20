%% Reader cylindrical coordinate 2D

clear all
ax = importdata('axis.txt');
%A0 = importdata('out0.txt');
A1 = importdata('out1.txt');

%% 
nr = 520%A3(1);
nz = 680%A3(2);

Ntime = 1%A3(3);
snap  = 10%A3(4);%1;%

Nsnap = snap;%floor(Ntime/snap);

r = ax(1:nr);
z = ax(nr+1:nz+nr);

[R Z] = meshgrid(r,z);

r0=0;


%% 

xmin = 0;
xmax = 24;
ymin =-15;
ymax = 15;

for j=1:Nsnap
    %clf
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
        'Color','w');
    
    
    PHI=reshape(A1(1+nr*nz*(j-1):nr*nz*j),nz,nr);
    
    %for jmovie=1:3
    % subplot(5,3,j)
    surfc(R,Z,log10(PHI+1e-12),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    ylabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    xlabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    hold on
    %plot3([0 max(r)],[0 0],[1 1]*(1),'Color','w')
    %plot3([r0 r0],[min(z) max(z)],[1 1]*(1),'Color','w')    
    
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String','| \phi |  ',...
     'fontsize',16)%,'fontweight','b');    
    %caxis([-12 0])
    title(['Time: ',num2str(1.*j),' a.u.'],'fontsize',16)
   % title(['Time: ',num2str(A2(j)),' a.u.'],'fontsize',16)
   axis tight
   % xlim([xmin xmax])
   % ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)
    
    %caxis([-6 -1])
  %  F = getframe(fig);
  %  aviobj = addframe(aviobj,F); 
%end 
    pause(0.2)
 %   display(j);
end    
    
%aviobj = close(aviobj);

%%
