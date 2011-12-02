clear all
A=importdata('out1.txt');
AU=importdata('out1U.txt');
r=importdata('ejes.txt');

rH=r(:,2);
rU=r(:,1);


B=reshape(A,1000,400);
BU=reshape(AU,1000,400);


figure
surf(BU,'FaceColor','interp','EdgeColor','none')
 camlight left; lighting phong
 axis tight
 %axis([1 500 1 500 -1 1])
 %caxis([-0.7 0.7])
 
 
 
 colorbar('location','south')
 view(2)
 
%%
figure
for i=1:100
    subplot(2,1,1)
    
    hold off
    plot(rU,log10(BU(:,i)),'LineWidth',3)
    hold on
    plot(rH,log10(B(:,i)),'r','LineWidth',3)
    axis tight
  %  xlim([400 600])
    
    subplot(2,1,2)
    
    phiinter=interp1(rU,BU(:,i),rH);

    phiinter2=interp1(rH,B(:,i),rU);
    
    hold off
    plot(rU,phiinter-BU(:,i),'g','LineWidth',3)
    hold on
    plot(r,phiinter2-B(:,i),'b','LineWidth',3)
    axis tight
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
 
 