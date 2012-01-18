%%%%%%  Pole Reader

A0=importdata('out0.txt');
A1=importdata('out1.txt');



nr=520;
nz=680;

Ntime=200;
Nsnap=30;

cc=Ntime/Nsnap;
%%

In=reshape(A0,nz,nr);
%Pole=reshape(P,nz,nr);

figure

surf(In,'FaceColor','interp','EdgeColor','none')
view(2)
axis tight

%%
% for i=1:1
%     Out=reshape(A1((i-1)*nz*nr+1:i*nz*nr),nz,nr);
%     
%     surf(Out,'FaceColor','interp',...
%         'EdgeColor','none')
%    
%     pause(0.3)
% end