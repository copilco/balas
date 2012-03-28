%% Reader 

clear all

A0 = importdata('out0.txt');
A1 = importdata('out1.txt');
A2 = importdata('out2.txt');
A3 = importdata('out3.txt');
A4 = importdata('out4.txt');
A5 = importdata('out5.txt');
A6 = importdata('out6.txt');

Aa = importdata('outa.txt');

%% Definition Grid
i = sqrt(-1);

x=A0;
dx=x(2)-x(1);
Nx=length(A0);

k  = A5;
dk = k(2)-k(1);
Nk=length(k);


Nsnap=A6(1,1);
NNx=length(A2)/Nsnap;
[X Index] = meshgrid(1:Nsnap,x);





%% Reorder data

clear PS
for j=1:Nsnap
    MSpectrum(:,j)= A4(1+(j-1)*Nk:j*Nk,3);%+i*A4(1+(j-1)*nk:j*nk,2);
    Asymmetry(:,j)= (rot90(MSpectrum(1:Nk/2,j),2)-MSpectrum(Nk/2+1:Nk,j))./(rot90(MSpectrum(1:Nk/2,j),2)+MSpectrum(Nk/2+1:Nk,j));
    
    AsymmetryDefTwo(j) = (dk*sum(rot90(MSpectrum(1:Nk/2,j).*(k(1:Nk/2).^2),2))-sum(MSpectrum(Nk/2+1:Nk,j).*(k(1+Nk/2:Nk).^2))*dk)...
        ./(sum(rot90(MSpectrum(1:Nk/2,j).*(k(1:Nk/2).^2),2))*dk+sum(MSpectrum(Nk/2+1:Nk,j).*(k(1+Nk/2:Nk).^2))*dk);
    PS(:,j) = A2(1+(j-1)*Nx:j*Nx);
    
end




%%
figure
plot(k(Nk/2+1:Nk),Asymmetry(:,Nsnap))
figure
plot(AsymmetryDefTwo)

%plot(k,MSpectrum(:,Nsnap-1))
%%

figure

for j=1:Nsnap
    clf
    plot(x(1:Nx),log10(PS(:,j)+1e-12),'r','linewidth',1)
    xlim([-200 200])
    pause(0.5)
end

%%


figure
ymin =-300;
ymax = 300;

surf(Index,X,log10(PS+1e-12),'FaceColor','interp','EdgeColor','none'); 
xlabel(' N/snap ','fontsize',16,...
        'fontname','Helvetica')%,'fontweight','b');
ylabel(' z (a.u.) ','fontsize',16,...
        'fontname','Helvetica')%,'fontweight','b');
    

    
h = colorbar('location','EastOutside');
set(get(h,'YLabel'),'String','log10(Prob)  ',...
    'fontsize',16)


view(2)
axis tight

xlim([ymin ymax])
h=gca;    
set(h,'fontsize',16)


%%
