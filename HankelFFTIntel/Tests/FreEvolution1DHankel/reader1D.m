A=importdata('outPhiRetrieved.txt');
r=importdata('rhofile.txt');
t=importdata('timefile.txt');

[T R]=meshgrid(t,r);

B=reshape(A,3000,20);
figure
surfc(T',R',log10(B'+1e-7),'FaceColor','interp',...
        'EdgeColor','none')
    view(2)
    
    grid on
    caxis([-7 1])
 %camlight left   
 colorbar
    axis tight
