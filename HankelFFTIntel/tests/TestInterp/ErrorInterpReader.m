%% Reader cylindrical coordinate 2D

clear all
ax = importdata('axis.txt');
%A0 = importdata('out0.txt');
%A1 = importdata('outH2U.txt');
%A2 = importdata('outU2H.txt');
A3 = importdata('outErrorH2U.txt');
A4 = importdata('outErrorU2H.txt');
%A5 = importdata('outDiffHankel.txt');
%A6 = importdata('outDiffCrank.txt');

%%
nr = 250;
nz = 500;

time=A3(:,1);

ErrorH=A3(:,2);
ErrorUin=A3(:,3);

ErrorU=A4(:,2);
ErrorHin=A4(:,3);

%% Plot the error in normalization H2U

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');

plot(time,log10(abs(ErrorH)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorUin)),'ro','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Comparision in error normalization. H2U','fontsize',16)
grid on
hleg1 = legend('Hankel','Uniform Interp');


%% Plot the error in normalization H2U

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],'Color','w');

plot(time,log10(abs(ErrorU)),'LineWidth',3)
hold on
plot(time,log10(abs(ErrorHin)),'ro','LineWidth',1)
hold off
xlabel('Temporal iterations','fontsize',12)
ylabel('Logaritmic error','fontsize',12)
title('Comparision in error normalization. U2H','fontsize',16)
grid on
hleg1 = legend('Uniform','Hankel Interp');

%% Plot the difference between points in logaritmic scale (2D density plot)

% Axes
r = ax(1:nr);
z = ax(nr+1:nz+nr);

[R Z] = meshgrid(r,z);

snap = length(A5)/nr/nz;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.5],'Color','w');

for k=1:snap
    
   
    subplot(2,2,[1 3])
    DiffHankel = reshape(A5((k-1)*nr*nz+1:k*nr*nz),nz,nr);
    surf(R,Z,log(abs(DiffHankel)),'FaceColor','interp','EdgeColor','none')
    view(2)
    axis tight
    %daspect([1,1,1])
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String','Logaritmic difference betweeen Hankel and interpoled one.',...
        'fontsize',16)
    xlabel('\rho axis','fontsize',16);
    ylabel('z axis','fontsize',16);
    
    subplot(2,2,[2 4])
    DiffCrank = reshape(A6((k-1)*nr*nz+1:k*nr*nz),nz,nr);
    surf(R,Z,log(abs(DiffCrank)),'FaceColor','interp','EdgeColor','none')
    view(2)
    axis tight
    %daspect([1,1,1])
    h = colorbar('location','EastOutside');
    set(get(h,'YLabel'),'String','Logaritmic difference between Uniform and interpolated one',...
        'fontsize',16);
    xlabel('\rho axis','fontsize',16);
    ylabel('z axis','fontsize',16);
    
    %title('Comparision')
    

end