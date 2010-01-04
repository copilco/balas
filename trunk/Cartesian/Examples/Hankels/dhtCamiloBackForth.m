gauss1D = @(x, x0, FWHM) exp(-2*log(2)*((x-x0)/FWHM).^2);
gauss1 = @(x) 2*x.*exp(-x.^2);

% Nr = 1024;			%	Number of sample points
% r_max = .05;		%	Maximum radius (5cm)
% dr = r_max/(Nr-1);	%	Radial spacing
% nr = (0:Nr-1)';		%	Radial pixels
% r = nr*dr;			%	Radial positions
% Dr = 5e-3;			%	Beam radius (5mm)



Rmax=5%2.25
Kmax=10%11
[H,k,r,I,K,R,h]=dht2(gauss1,Rmax)
%[H,k,r,I,h]=fht2(gauss1,Rmax,Kmax)

A=importdata('out0.txt');
B=importdata('out1.txt');
C=importdata('out2.txt');


%figure
subplot(2,3,4)
%plot(sqrt(A(:,1).^2+A(:,2).^2))
plot(A(:,1),A(:,2))
title('kernel22')

subplot(2,3,5)
%plot(sqrt(B(:,1).^2+B(:,2).^2))
plot(B(:,1),B(:,2))
title('tempC++')

subplot(2,3,6)
%plot(C(:,1),C(:,2))
plot(fftshift(sqrt(C(:,1).^2+C(:,2).^2)))
title('tempC++2')


D=importdata('out3.txt');
E=importdata('out4.txt');
F=importdata('out5.txt');

figure
subplot(1,3,1)
plot(D(:,1),D(:,2),'LineWidth',4)
hold on
plot(r,h,'g')
title('r vs h')
 

 subplot(1,3,2)
% NN=length(h)
 plot(E(:,1),E(:,2),'LineWidth',4)
% %plot(B(:,2))
 hold on
 plot(k,H,'g')
 title('k vs H')


 


%[aa,II,horig]=ifht2(H,k,r)
%figure
subplot(1,3,3)
plot(F(:,1),F(:,2)/pi/pi,'LineWidth',4)
hold on
%plot(r,horig,'g')
title('r vs h')


 



% 
% subplot(3,1,1)
% plot(r,gauss1D(r, 0, 1))
% xlim([0 Rmax])
% 
% subplot(3,1,2,'r')
% plot(k,H)
% 
% subplot(3,1,1)
% hh=ifht(H,k,r,0);
% hold on
% plot(r,hh,'g')


