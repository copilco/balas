%%%%%%%%%   New Hankel example

clear all

%% Set parameters

Nr = 100;			%	Number of sample points
r_max = 200;       %	Maximum radius (5cm)
dt=0.5;
%dr = r_max/(Nr-1);	%	Radial spacing
%nr = (0:Nr-1)';		%	Radial pixels
%r = nr*dr;			%	Radial positions
%Dr = 5e-3;			%	Beam radius (5mm)
%Kr = 5000;			%	Propagation direction
%Nz = 200;			%	Number of z positions
%z_max = .25;		%	Maximum propagation distance
%dz = z_max/(Nz-1);
%z = (0:Nz-1)'*dz;	%	Propagation axis
Nloops=100;

time=1:Nloops;

%% Load HankelMatrix

H = hankel_matrix(0,r_max,Nr);

K = 2*pi*H.V;

%% Set the differential
ddr=zeros(Nr,1);

for i=1:Nr-1
    ddr(i)=H.r(i+1)-H.r(i);
end
ddr(Nr)=r_max-H.r(Nr);



ddv=zeros(Nr,1);

for i=1:Nr-1
    ddv(i)=H.v(i+1)-H.v(i);
end
ddv(Nr)=H.V-H.v(Nr);

%figure
%plot(ddr,'o')

%% Plot the initial wavefuntion
 % Choose one sigma
%sigma=.7;
sigma=10;

phir=exp(-(H.r-r_max/2.).*(H.r-r_max/2.)/sigma/sigma);%(r_max/500)/(r_max/500));


norm=0.;
for i=1:Nr
    norm=norm+ddr(i)*H.r(i)*real(conj(phir(i))*phir(i));
end


phir=phir/sqrt(norm);


norm=0.;
for i=1:Nr
    norm=norm+ddr(i)*H.r(i)*real(conj(phir(i))*phir(i));
end

disp('Original norm:')
disp(norm)

figure
plot(H.r,abs(phir),'-o')
ylim([0 .1])
xlabel('r')
ylabel('Función de onda')
title('Función de onda inicial')

%% Plot the initial spectrum
phik = qdht(phir./H.JR , H, 128).*H.JV;

scrsz=get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
plot(H.v,abs(phik),'b')
hold on
plot(H.v,phik,'r-o')
hold off
%xlim([0 .5])
xlabel('Frequencies')
ylabel('Intensity (arbitrary units)')
title('Wavefunction in frequency space')

norm2=0.;
    for k=1:Nr
        norm2=norm2+ddv(k)*H.v(k)*real(conj(phik(k))*phik(k));
    end
    
disp('Original norm in frequency space:')
disp(norm2);

%% Wavefuntion evolution
total_norm=zeros(Nloops,1);
C=zeros(Nr,Nloops);
disp('Performing Hankel transform ...');
scrsz=get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])


for i=1:Nloops
    phik = qdht(phir./H.JR , H, 128 );
    
    fase=4.*pi*pi*H.v.*H.v.*dt/2.;
    phik=phik.*exp(j*fase);
    
    PHIRZ = iqdht(phik, H,128);
    
    phir=PHIRZ.*H.JR;
    
    C(:,i)=phir;
    
    plot(abs(phir),'r' )
    ylim([0 .1])
    xlabel('r')
    ylabel('Función de onda')
    title('Función de onda propagada')

    norm1=0.;
    for k=1:Nr
        norm1=norm1+ddr(k)*H.r(k)*real(conj(phir(k))*phir(k));
    end
    
    
    disp('Retreived error: ')
    disp(1.-norm1)
    
    total_norm(i,1)=norm1;
    
    pause(0.1)
    
end

%% Spectrum evolution
total_norm2=zeros(Nloops,1);
scrsz=get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.2 scrsz(4)/1.2])

for i=1:Nloops
    PHIK = qdht(phir./H.JR , H, 128 );
    
    fase=2.*pi*2.*pi*H.v.*H.v*dt/2.;
    PHIK=PHIK.*exp(j*fase);
    
    phik=PHIK.*H.JV;
    
    PHIRZ = iqdht(PHIK, H,128);
    
    phir=PHIRZ.*H.JR;
    
    
    
    plot(H.v,(phik),'r-o')
    hold on
    plot(H.v,abs(phik),'b')
    hold off
    %xlim([0 1])
    ylim([-50 50])
    xlabel('v')
    ylabel('Wavefuntion')
    title('Propagated wavefunction in frequency space')

    norm2=0.;
    for k=1:Nr
        norm2=norm2+ddv(k)*H.v(k)*real(conj(phik(k))*phik(k));
    end
     
    
    disp('Retreived error: ')
    disp(1.-norm2)
    
    total_norm2(i,1)=norm2;
    
    pause(0.1)
    
end

%% Final snapshot

figure
surf(time,H.r,abs(C),'FaceColor','interp','EdgeColor','none')
view(2)
axis tight
xlabel('Pasos temporales')
ylabel('Coordenada r')



%% Error in coordinate space

resta=1.-total_norm;
figure
subplot(2,1,1)
plot(resta,'o')
xlabel('Pasos temporales')
title('Error en la norma en coordenadas')
subplot(2,1,2)
plot(abs(log10(resta)))
xlabel('Pasos temporales')
title('Log10 norm error')

%% Error in frequency space

resta2=1.-total_norm2;
figure
subplot(2,1,1)
plot(resta2,'o')
xlabel('Pasos temporales')
title('Error en la norma en frecuencias')
subplot(2,1,2)
plot(abs(log10(resta2)))
xlabel('Pasos temporales')
title('Log10 norm error')

