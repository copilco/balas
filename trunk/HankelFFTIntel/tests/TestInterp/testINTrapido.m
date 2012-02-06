%%%%%%  Pole Reader

H=importdata('outH2U.txt');

nr=500;
nz=750;

Ntime=100;
Nsnap=20;
snap=Ntime/Nsnap;
%%



for k=1:15%snap
    
    Hank=reshape(H((k-1)*nr*nz+1:k*nr*nz,2),nz,nr);
    Interp=reshape(H((k-1)*nr*nz+1:k*nr*nz,3),nz,nr);
    
    
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
   
    subplot(1,2,1)
    surf(Hank,'FaceColor','interp','EdgeColor','none')
    view(2)
    axis tight

    subplot(1,2,2)
    %scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    surf(Interp,'FaceColor','interp','EdgeColor','none')
    view(2)
    axis tight

end