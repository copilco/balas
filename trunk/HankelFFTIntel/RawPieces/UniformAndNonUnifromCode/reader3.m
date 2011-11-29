%% Reader Non Uniform Cylindrical Crank-Nicholson 1D
% By Alexis Chacon November, 2011

clear all
A0=importdata('out0.txt');
A1=importdata('out1.txt');
A2=importdata('out2.txt');
A3=importdata('out3.txt');



%% Parameters and axes
snap  = length(A2)-2;
Ntime = A2(snap+1);
snapp = floor(Ntime/snap);
dt    = A2(snap+2); 
t     = A2(1:snap);


Nr    = length(A0(:,1));
dr    = A0(2,1)-A0(1,1);
r     = A0(:,1);

rw    = 15;
r0    = 185;
r1    = 100;

r_n0  = floor( (r1 - rw/2)/dr )+1;
r_n1  = floor( (r1 + rw/2)/dr )+1;

DensProb = reshape(A1,Nr,snap);



[T R] = meshgrid(t,r);


%% Ploting graph

xmin = 0;
xmax = 400;
ymax=2.4;
ymin=0;
 

for ktime=1:snap
    clf
    
for jmovie=1:1    
    aa=area(r(1:Nr),DensProb(1:Nr,ktime)+0.);    
    set(aa(1),'FaceColor',[0.1 0.1 1])
    set(aa(1),'EdgeColor',[0.1 0.1 1])
    
    xlabel('\rho (a.u.) ','fontsize',16 )
    ylabel('| \psi(\rho) | ','fontsize',16 )
    title(...%['norm: '  ,num2str(norma(r,DensProb(:,ktime))),...
        ['time: '  ,num2str(t(ktime)),' a.u.'...
        ';  ktime: ' ,num2str(ktime*snapp)],'fontsize',16)    
    
    hold on

%     aa=area(r(r_n0+1:r_n1),DensProb(r_n0+1:r_n1,ktime)+0.);    
%     set(aa(1),'FaceColor',[0.4 0.4 0.4])
%     set(aa(1),'EdgeColor',[0.4 0.4 0.4])
% 
%     aa=area(r(r_n1+1:Nr),DensProb(r_n1+1:Nr,ktime)+0.);    
%     set(aa(1),'FaceColor',[0.1 0.9 0])
%     set(aa(1),'EdgeColor',[0.1 0.9 0])    
    
    fac=1;%max(A0(:,2))*max(max(DensProb));
    plot(r,A0(:,2)*fac+0,'r','linewidth',1)
    if (mod(ktime,3)==0)
       line([r1 r1],[-2 2],...
           'Color','r','linewidth',2,'linestyle','--') 
       %line([192.4+15 192.4+15],[-2 2],...
        %   'Color','g','linewidth',2,'linestyle','--')        
    end
    ylim([ymin ymax])
    xlim([xmin xmax])
    
    h=gca;    
    set(h,'fontsize',16)
    %F = getframe(fig);
    %aviobj = addframe(aviobj,F);
end
    pause(0.01)
end

%aviobj = close(aviobj);



%% Spacing vs numer of point in rho; 
%  and rho vs number of point in rho

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.7 scrsz(4)*0.5],...
         'Color','w');

subplot(2,1,1)
plot(A0(:,3))
ylabel('\delta\rho (a.u. ) ','fontsize',16 )
xlabel('pixel ','fontsize',16 )
h=gca;    
set(h,'fontsize',16)
axis tight

subplot(2,1,2)
plot(r,'r')
ylabel('\rho (a.u. ) ','fontsize',16 )
xlabel('pixel ','fontsize',16 )
h=gca;    
set(h,'fontsize',16)
axis tight



%% Surfc plot of the wavefunction progation

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)*0.9 scrsz(4)*0.5],...
    'Color','w');

ymax = 200;
xmax = max(t);
xmin = 0.0;

surfc(T,R,log10(DensProb+1e-7),...
        'FaceColor','interp',...
        'EdgeColor','none')    
view(2)
axis tight
xlabel('Time ','fontsize',16 )
ylabel('\rho (a.u.) ','fontsize',16 )
title('NonUniformGrid drj = 1*10^{-5}*j + dr; dr=0.1 ','fontsize',16 )

colorbar('location','EastOutside');
set(gca,'fontsize',16)

%    set(get(h,'XLabel'),'String','Prob  ',...
%     'fontsize',14,'fontweight','b');
%xlim([min(delay) max(delay)])

h=gca;    
set(h,'fontsize',16)

ylim([0 ymax])    
xlim([xmin xmax]) 
caxis([-7 1])


%% Numeric Error on Norm vs Time Iteration
figure
plot((1:snap)*snapp,log10(abs(A3-1)+1e-17))
ylabel('RelativeErrorNorm ','fontsize',16 )
xlabel('Time Iteration ','fontsize',16 )
title('NonUniformGrid drj = 1*10^{-5}*j + dr; dr=0.1 ','fontsize',16 )
h=gca;    
set(h,'fontsize',16)
%%
