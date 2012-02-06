%% Reader cylindrical coordinate 2D, Imaginary time propagation

clear all

A0 = importdata('out0.txt');
A1 = importdata('out1.txt');
%A2 = importdata('out2.txt');
A3 = importdata('out3.txt');
A4 = importdata('out4.txt'); 
A5 = importdata('out5.txt'); 


%% Grid parameters
nr = 250; %A3(1);
nz = 400; %A3(2);


snap  = length(A3)/nr/nz;%1;%
Nsnap = snap-1;%floor(Ntime/snap);

r = A0(1:nr);
z = A0(nr+1:nz+nr);

[Z R] = meshgrid(z,r);

dt = 0.05;

t     = A5(:,1);
Exuv  = A5(:,2);
ksnap = floor(length(t)/snap);



%% Movie Convergence bound state



xmin = -25.;
xmax =  25.;
ymin =  0. ;
ymax =  20.;


for j=1:10:Nsnap
    
    %scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.8],...
    %    'Color','w');    
    figure
    
    PHI  = reshape(A3(1+nr*nz*(j-1):nr*nz*j),nz,nr);    
    
    
    surfc(Z,R,log10(PHI'+1e-12),...
        'FaceColor','interp',...
        'EdgeColor','none')
    
    view(2)
    axis tight
    
    h=gca;    
    set(h,'fontsize',16)
    
    xlabel('z (a.u.) ','fontsize',16)%,'fontweight','b')
    ylabel('\rho (a.u.) ','fontsize',16)%,'fontweight','b')  
    
    
    
    
    title(j)
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String',' |\rho \phi |^2  ',...
        'fontsize',16)%,'fontweight','b');    
    
    caxis([-12 0])
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    h=gca;    
    set(h,'fontsize',16)
    
    
    
    pause(0.01)

end    
    


%% Convergence error

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.5],...
        'Color','w');

semilogy((1:40:4000),A4(:,2),'linewidth',3)
xlabel('iteration time ','fontsize',16)%,'fontweight','b')
ylabel(' energy error  ','fontsize',16)

set(gca,'fontsize',16)

%%



