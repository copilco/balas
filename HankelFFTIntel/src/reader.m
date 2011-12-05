clear all
A=importdata('out1.txt');
AU=importdata('out1U.txt');
r=importdata('ejes.txt');

rH=r(:,2);
rU=r(:,1);


B=reshape(A,1000,100);
BU=reshape(AU,1000,100);


figure
surf(BU,'FaceColor','interp','EdgeColor','none')
camlight left; lighting phong
axis tight
%axis([1 500 1 500 -1 1])
%caxis([-0.7 0.7])
view(2)

figure
surf(B,'FaceColor','interp','EdgeColor','none')
camlight left; lighting phong
axis tight
%axis([1 500 1 500 -1 1])
%caxis([-0.7 0.7])
view(2)



colorbar('location','south')
view(2)

%%
figure
for i=1:100
    
    
    
    subplot(3,1,1)
    
    hold off
    plot(rU,log10(BU(:,i)),'LineWidth',3)
    hold on
    plot(rH,log10(B(:,i)),'r','LineWidth',3)
    axis tight
    xlabel('RU (a.u.), RH (a.u.)')
    ylabel('Log dP')
    %  xlim([400 600])
    
    h = legend('Uniform','Hankel',2);
    set(h,'Interpreter','none','Location','NorthEast')
    
    set(h,'fontsize',14,'FontWeight','b');
     set(gca,'fontsize',14,'FontWeight','b');
    
    
    subplot(3,1,2)
    
    hold off
    plot(rU,(BU(:,i)),'LineWidth',3)
    hold on
    plot(rH,(B(:,i)),'r','LineWidth',3)
    axis tight
    %  xlim([400 600])
    
    xlabel('RU (a.u.), RH (a.u.)')
    ylabel('Log dP')
    
    h = legend('Uniform','Hankel',2);
    set(h,'Interpreter','none','Location','NorthEast')
    
    
    set(h,'fontsize',14,'FontWeight','b');
    
    set(gca,'fontsize',14,'FontWeight','b');
    
    
    subplot(3,1,3)
    
    phiinter=interp1(rU,BU(:,i),rH);
    
    phiinter2=interp1(rH,B(:,i),rU);
    
    hold off
    plot(rU,phiinter-BU(:,i),'g','LineWidth',3)
    hold on
    plot(r,phiinter2-B(:,i),'b','LineWidth',3)
    axis tight
    
     h = legend('Uniform','Hankel',2);
    set(h,'Interpreter','none','Location','NorthEast')
    
    xlabel('RU (a.u.), RH (a.u.)')
    ylabel('dPU - dPH')
    
    set(h,'fontsize',14,'FontWeight','b');
     set(gca,'fontsize',14,'FontWeight','b');
    pause(0.3)
end


figure
surf(B-BU,'FaceColor','interp','EdgeColor','none')
camlight left; lighting phong
axis tight
%axis([1 500 1 500 -1 1])
%caxis([-0.7 0.7])



colorbar('location','south')
view(2)


figure
plot(rU,rH)
