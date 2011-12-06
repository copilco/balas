% A=importdata('outPhiRetrieved.txt');
% r=importdata('rhofile.txt');
% t=importdata('timefile.txt');
% 
% [T R]=meshgrid(t,r);
% 
% B=reshape(A,10000,100);
% figure
% surfc(T',R',log10(B'+1e-7),'FaceColor','interp',...
%         'EdgeColor','none')
%     view(2)
%     
%     grid on
%     caxis([-7 1])
%  %camlight left   
%  colorbar
%     axis tight
A=importdata('outPhiRetrieved.txt');



B=reshape(A,1000,100);
figure
surf(B,'FaceColor','interp','EdgeColor','none')
 camlight left; lighting phong
 axis tight
 %axis([1 500 1 500 -1 1])
 %caxis([-0.7 0.7])
 
 
 
 colorbar('location','south')
 view(2)