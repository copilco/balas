clear all

Nr=1000
fid = fopen('wUbalas.bin','r');
AU=fread(fid,Nr,'*double');

fidH = fopen('wHbalas.bin','r');
AH=fread(fidH,Nr,'*double');

fid1= fopen('EjeUbalas.bin','r');
rU=fread(fid1,Nr,'*double');

fid2= fopen('EjeHbalas.bin','r');
rH=fread(fid2,Nr,'*double');


%fstream laphi ("wUbalas.bin", ios::out | ios::binary);
%	laphi.read ((char*)&temp[0], sizeof (double)*(w.Nr) );
%	laphi.close();
	
%	fstream raxis ("EjeUbalas.bin", ios::out | ios::binary);
%	raxis.read ((char*)&temp[0], sizeof (double)*(w.Nr) );
%	raxis.close();
    
    
%A=importdata('out1.txt');
%AU=importdata('out1U.txt');
%r=importdata('ejes.txt');

%rH=r(:,2);
%rU=r(:,1);



%%
figure
       
    subplot(3,1,1)
    
    hold off
    plot(rU,log10(AU),'LineWidth',3)
    hold on
    plot(rH,log10(AH),'r','LineWidth',3)
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
    plot(rU,AU,'LineWidth',3)
    hold on
    plot(rH,AH,'r','LineWidth',3)
    axis tight
    %  xlim([400 600])
    
    xlabel('RU (a.u.), RH (a.u.)')
    ylabel('Log dP')
    
    h = legend('Uniform','Hankel',2);
    set(h,'Interpreter','none','Location','NorthEast')
    
    
    set(h,'fontsize',14,'FontWeight','b');
    
    set(gca,'fontsize',14,'FontWeight','b');
    
    
    subplot(3,1,3)
    
    phiHMat=interp1(rU,AU,rH,'spline');
    phiUMat=interp1(rH,AH,rU,'spline');
    
  %  phiinter2=interp1(rH,A,rU);
    
    hold off
    plot(rU,phiUMat-AU,'b*-','LineWidth',3)
    hold on
    plot(rH,phiHMat-AH,'go','LineWidth',3)
    axis tight
    %xlim([100 110])
     h = legend('Uniform','Hankel',2);
    set(h,'Interpreter','none','Location','NorthEast')
    
    xlabel('RU (a.u.), RH (a.u.)')
    %ylabel('dPU - dPH')
    
    set(h,'fontsize',14,'FontWeight','b');
     set(gca,'fontsize',14,'FontWeight','b');
    pause(0.3)

    
    
    fid = fopen('phiHankelmatlab.bin', 'w');
    fwrite(fid, phiHMat, 'double');
    fclose(fid)
    
    fidH = fopen('phiUniformmatlab.bin', 'w');
    fwrite(fidH, phiUMat, 'double');
    fclose(fidH)
    
%     %     phiH=interp1(rU,AU,rH);
% %     
%      aa=isnan(phiHMat);
% %     
%      for i=1:Nr
%         if(aa(i)==1)     
%          phiHMat(i)=0;
%          i
%          end
%      end
%      figure
%       aa=isnan(phiHMat);
%      plot(aa)
% %     
%      fid = fopen('phiHankelmatlab.bin', 'w');
%      fwrite(fid, phiHMat, 'double');
%      fclose(fid)
%     
%      
%        aa=isnan(phiHMat);
% %     
%      for i=1:Nr
%         if(aa(i)==1)     
%          phiH(i)=0;
%          i
%          end
%      end
%      figure
%       aa=isnan(phiHMat);
%      plot(aa)
% %     
%      fid = fopen('phiHankelmatlab.bin', 'w');
%      fwrite(fid, phiHMat, 'double');
%      fclose(fid)
%      
%     %%
%     
%     
%     
%     
% 
