%% Reader 

clear all

A1=importdata('out1.txt');
%A1=importdata('out2.txt');


%% 
i=sqrt(-1);
figure
subplot(2,1,1)
plot(A1(:,1),unwrap(angle(A1(:,2)+i*A1(:,3))))
subplot(2,1,2)
plot(A1(:,1),A1(:,4),'r')

%%
figure
nk=length(A1(:,1));
plot((rot90(A1(1:nk/2,2),2)-A1(nk/2+1:nk,2))./(rot90(A1(1:nk/2,2),2)+A1(nk/2+1:nk,2)))
set(gca,'fontsize',16)

%%