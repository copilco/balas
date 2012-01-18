%%Comparing Hankel matrix
clear all
R=200.
Nr=1000



 fid = fopen('BesselMatrix.bin');
 MH = fread(fid,[Nr Nr],'double');
 fclose(fid);


MHMat = hankel_matrix(0,R,Nr)


resta=MH-MHMat.T;% Matmat.T-Creal;


 figure

 subplot(2,2,1)
 surf(log10(abs(resta)),'FaceColor','interp','EdgeColor','none')
 title('Error log10(abs(HankeMat-HankelMat.T))')

view([-208 78])


%%

fid = fopen('HankelRhoaxis.bin');
rho = fread(fid,Nr,'double');
fclose(fid);


resta=rho-MHMat.r;% Matmat.T-Creal;

subplot(4,4,[3 4])
plot(rho,'g')
hold on
plot(MHMat.r,'bo','MarkerSize',4 )
title('rho Axis (green) v(blue) ))')


subplot(4,4,[7 8])
plot(log10(abs(resta)),'g')




subplot(4,4,[15 16])
fid = fopen('HankelVaxis.bin');
v = fread(fid,Nr,'double');
fclose(fid);

resta=v-MHMat.v;% Matmat.T-Creal;

hold on
plot(log10(abs(resta)),'b')
title('Error rho (green) v(blue) ))')

subplot(4,4,[11 12])
plot(v,'g')
hold on
plot(MHMat.v,'bo','MarkerSize',4 )
title('V Axis (green) v(blue) ))')


 
 %%
 
 
 

fid = fopen('Hankelm1.bin');
m1 = fread(fid,Nr,'double');
fclose(fid);


resta=m1-MHMat.JR;% Matmat.T-Creal;


fid2 = fopen('Hankelm2.bin');
m2 = fread(fid2,Nr,'double');
fclose(fid2);


resta2=m2-MHMat.JV;% Matmat.T-Creal;


 subplot(2,2,3)
 hold on
 plot(log10(abs(resta)),'g')
 
  plot(log10(abs(resta2)),'b')
 title('Errot m1 (green) m2 (blue)')
 
 
 
 

%resta=MH-MHMat.JV;% Matmat.T-Creal;
% %%
% 
% figure
% subplot(2,1,1)
% plot(resta)
% title('resta eje m1')
% subplot(2,1,2)
% plot(log10(abs(resta)),'r')
% title('log10(abs(resta))')
% %subplot(2,2,4)
% %surf(abs(log10(resta)))
% %title('abs(log10(resta))')
% 

%%

% figure
% subplot(2,1,1)
% surf(resta,'FaceColor','interp','EdgeColor','none')
% title('resta')
% subplot(2,2,3)
% surf(log10(abs(resta)),'FaceColor','interp','EdgeColor','none')
% title('log10(abs(resta))')
% subplot(2,2,4)
% surf(abs(log10(resta)),'FaceColor','interp','EdgeColor','none')
% title('abs(log10(resta))')

%%
% 
% figure
% plot(MHMat.r,'b')
% hold on
% plot(raxis,'ro')
% 
% %%
% 
% restar=Matmat.r-raxis;
% 
% figure
% subplot(2,1,1)
% plot(restar)
% title('Resta del eje r')
% subplot(2,1,2)
% plot(log10(abs(restar)))
% title('Log10 de la resta del eje r')