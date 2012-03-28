%% Reader 

clear all

A0=importdata('out3.txt');
A1=importdata('out2.txt');


%% 

nx=length(A0);
figure
Nsnaphop=length(A1)/nx;
dx=A0(2,1)-A0(1,1)

for j=1:Nsnaphop
   clf
   plot(A0(:,1),A0(:,2),A0(:,1),A1(1+(j-1)*nx:j*nx)/dx-0.5)
   title(j)
   xlim([-50 50])
   pause(1)
end

%%
for j=1:Nsnaphop

    normA(j)=sum(A1(1+(j-1)*nx:j*nx));

end

%%

figure
plot(log10(normA))
set(gca,'fontsize',16)
%%

%%