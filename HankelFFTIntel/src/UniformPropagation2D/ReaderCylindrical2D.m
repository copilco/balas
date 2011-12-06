%% Reader cylindrical coordinate 2D

clear all
A0 = importdata('out0.txt');
A1 = importdata('out2.txt');
%A2 = importdata('out2.txt');
%A3 = importdata('out3.txt');

%% 
nr = 100%A3(1);
nz = 200%A3(2);

Ntime = 1%A3(3);
snap  = 12%A3(4);%1;%

Nsnap = snap;%floor(Ntime/snap);

r = A0(1:nr);
z = A0(nr+1:nz+nr);

[R Z] = meshgrid(r,z);

r0=0;


%% 

%aviobj = avifile('GEWP_C2D.avi');
scrsz = get(0,'ScreenSize');
fig=figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.5],...
    'Color','w');
xmin = 0;
xmax = 24;
ymin =-15;
ymax = 15;
for j=1:Nsnap
   % clf
    figure
    PHI=reshape(A1(1+nr*nz*(j-1):nr*nz*j),nz,nr);
    
%for jmovie=1:3
   % subplot(5,3,j)
    surfc(R,Z,log10(PHI+1e-8),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    ylabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    xlabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    hold on
    plot3([0 max(r)],[0 0],[1 1]*(1),'Color','w')
    plot3([r0 r0],[min(z) max(z)],[1 1]*(1),'Color','w')    
    
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
