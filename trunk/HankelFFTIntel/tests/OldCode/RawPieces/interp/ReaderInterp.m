%% Reader
in0=importdata('in0.txt');
out0=importdata('out0.txt');
out1=importdata('out1.txt');

%%
figure
%subplot(2,1,1)
plot(in0(:,1),in0(:,2),'r-')
hold on
plot(out0(:,1),out0(:,2),'bo')


xlabel('r (a.u.)');
ylabel('\Psi (en unidades arbitrarias)');
title('Interpolación. En rojo la original')

%%
Nr=1000;
Nt=5000;
figure

for i=1:Nt
    outA=out1((i-1)*Nr+1:i*Nr,2);
    outB=out1((i-1)*Nr+1:i*Nr,3);

    plot(outA,outB)
    ylim([0 .1])
    %axis tight
    pause(0.01)
    
end  

%%

C=reshape(out1(:,3),Nr,Nt);
time=1:Nt;
pts=1:Nr;
%%

surfc(C,'FaceColor','interp','EdgeColor','none')
view(2)
axis tight
view(2)
axis tight
xlabel('Temporal step')
ylabel('Number of points')
title('Interpoled Crank Evolution')