%% Reader 

clear all

A0=importdata('out0.txt');
A1=importdata('out1.txt');


%% 

nx=length(A0);
figure
Nsnaphop=length(A1)/nx;
for j=1:Nsnaphop
   clf
   plot(A0,-1./sqrt(2.+A0.*A0),A0,A1(1+(j-1)*nx:j*nx).^1-0.5)
   title(j)
   xlim([-50 50])
   pause(0.1)
end

%%