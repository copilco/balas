%%%%%%  Pole Reader

A0=importdata('out0.txt');
A1=importdata('out1.txt');

nr=500;
nz=500;

Ntime=400;
Nsnap=30;

cc=Ntime/Nsnap;
%%

In=reshape(A0,nz,nr);
%Pole=reshape(P,nz,nr);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
    'Color','w');

surf(In,'FaceColor','interp','EdgeColor','none')
view(2)
axis tight

%%

% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
%     'Color','w');
% 
% 
% for i=1:20
%     
%     clf;
%     
%     
%     Out=reshape(A1((i-1)*nz*nr+1:i*nz*nr),nz,nr);
%     
%     surf(Out,'FaceColor','interp',...
%         'EdgeColor','none')
%     view(2)
%     axis tight
% 
%     pause(0.3)
% end