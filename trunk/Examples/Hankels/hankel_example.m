%%	Gaussian function
gauss1D = @(x, x0, FWHM) exp(-2*log(2)*((x-x0)/FWHM).^2);

%%	Initialise grid
disp('Initialising data ...');
Nr = 1024;			%	Number of sample points
r_max = .05;		%	Maximum radius (5cm)
dr = r_max/(Nr-1);	%	Radial spacing
nr = (0:Nr-1)';		%	Radial pixels
r = nr*dr;			%	Radial positions
Dr = 5e-3;			%	Beam radius (5mm)
Kr = 5000;			%	Propagation direction
Nz = 200;			%	Number of z positions
z_max = .25;		%	Maximum propagation distance
dz = z_max/(Nz-1);
z = (0:Nz-1)'*dz;	%	Popagation axis

%%	Setup Hankel transform structure
disp('Setting up Hankel transform structure ...');
H = hankel_matrix(0, r_max, Nr);
K = 2*pi*H.V;		%	Maximum K vector

%%	Generate electric field:
disp('Generating electric field ...');
Er = gauss1D(r, 0, Dr).*exp(i*Kr*r);	%	Initial field
ErH = spline(r, Er, H.r);				%	Resampled field 

%%	Perform Hankel Transform
disp('Performing Hankel transform ...');
EkrH = qdht(ErH, H, 3);		%	Convert from physical field to physical wavevector

%%	Propagate beam
disp('Propagating beam ...');
EkrH_ = EkrH./H.JV;		%	Convert to scaled form for faster transform
phiz = @(z) (sqrt(K^2 - H.kr.^2) - K)*z;	%	Propagation phase
EkrHz = @(z) EkrH_ .* exp(i*phiz(z));		%	Apply propagation
ErHz = @(z) iqdht(EkrHz(z), H);				%	iQDHT (no scaling)
Erz = @(z) spline(H.r, ErHz(z).*H.JR, r);	%	Interpolate & scale output

Irz = zeros(Nr, Nz+1);
Irz(:, 1) = abs(Er).^2;
for n=1:Nz
	Irz(:, n+1) =abs(ErHz(z(n))).^2; %% imag(ErHz(z(n))); %%abs(ErHz(z(n))).^2;
end

%%	Plot field
figure(1)
subplot(2,1,1)
plot(r*1e3, [abs(Er).^2 unwrap(angle(Er))], ...
	H.r*1e3, [abs(ErH).^2 unwrap(angle(ErH))], '+');
title('Initial electric field distribution')
xlabel('Radial co-ordinate (r) /mm');
ylabel('Field intensity /arb.');
legend('|E(r)|^2', '\phi(r)', '|E(H.r)|^2', '\phi(H.r)');
axis([0 10 0 1]);

subplot(2,1,2)
plot(H.kr, abs(EkrH).^2)
title('Radial wave-vector distribution');
xlabel('Radial wave-vector (k_r) /rad m^{-1}');
ylabel('Field intensity /arb.');
axis([0 1e4 0 max(abs(EkrH).^2)]);

figure(2)
imagesc(z, r*1e3, Irz);
colorbar;
set(gca, 'YDir', 'normal');
title('Radial field intensity as a function of propagation for annular beam');
xlabel('Propagation distance (z) /m');
ylabel('Radial position (r) /mm');
axis([0 .25 0 25])

disp('Complete!');

figure
surf(Irz(:,2:Nz+1),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
%daspect([5 5 1])
axis tight
%view(2)
%camlight left

AA=importdata('pout5.txt');
AAA=reshape(AA,Nr,Nz);
figure
surf(AAA,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
%daspect([5 5 1])
axis tight
%view(2)
%camlight left


% figure
% bb=importdata('pout3.txt');
% plot(bb(:,2),bb(:,3),'.r')
% hold on
% plot(EkrH_)


figure
a0=importdata('pout0.txt');
plot(a0(:,3),a0(:,4),'.r')
hold on
plot(ErH,'o','MarkerSize',10)
title('Input Function')

figure
a1=importdata('pout1.txt');
plot(a1(:,2),a1(:,3),'.r')
hold on
plot(EkrH)
title('f2 vs EKrH')

figure
a2=importdata('pout2.txt');
plot(a2(:,2),a2(:,3),'.r')
hold on
plot(EkrH_)
title('F2 vs EKrH_')

figure
a4=importdata('pout4.txt');
plot(a4(:,2),a4(:,3),'.r')
hold on
plot(ErHz(100*dz),'o','MarkerSize',10)
title('Fretrieved vs ErHz')



 figure
 bb=importdata('pout3.txt');
 plot(bb(:,2),bb(:,3),'.r')
 hold on
 plot(EkrH_ .* exp(i*phiz(100*dz)),'o','MarkerSize',10 )