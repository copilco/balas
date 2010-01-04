figure
A0=importdata('pout0.txt');
Nr=1024
Nt=120
B0=reshape(A0(:,2),Nt,Nr);

surfc(B0,'FaceColor','interp',...
          'EdgeColor','none')
 grid off
 view([0 90])
 colorbar
% xlabel('z (a.u.)','FontSize',12,'FontWeight' ,'bold')
% ylabel('Time (a.u.)','FontSize',12,'FontWeight','bold')
% caxis([0 0.4])
 axis tight
 
figure
B00=reshape(A0(:,2),Nr,Nt);
surfc(B00,'FaceColor','interp',...
          'EdgeColor','none')
 grid off
 view([0 90])
 colorbar
 xlabel('z (a.u.)','FontSize',12,'FontWeight' ,'bold')
% ylabel('Time (a.u.)','FontSize',12,'FontWeight','bold')
% caxis([0 0.4])
 axis tight
 

figure
A1=importdata('pout1.txt');
Nr=1024
Nt=120
B1=reshape(A1(:,1),Nt,Nr);

surfc(B1,'FaceColor','interp',...
          'EdgeColor','none')
 grid off
 view([0 90])
 colorbar
% xlabel('z (a.u.)','FontSize',12,'FontWeight' ,'bold')
% ylabel('Time (a.u.)','FontSize',12,'FontWeight','bold')
% caxis([0 0.4])
 axis tight
 
 figure
A2=importdata('pout2.txt');
Nr=1024
Nt=120
B2=reshape(A2(:,1),Nt,Nr);

surfc(B2,'FaceColor','interp',...
          'EdgeColor','none')
 grid off
 view([0 90])
 colorbar
% xlabel('z (a.u.)','FontSize',12,'FontWeight' ,'bold')
% ylabel('Time (a.u.)','FontSize',12,'FontWeight','bold')
% caxis([0 0.4])
axis tight