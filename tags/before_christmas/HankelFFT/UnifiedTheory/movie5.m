A=importdata('spout0.txt');


Nr=424
Nt=120



figure
for i=1:40
    B=A(i*Nr*Nt+1:(i+1)*Nr*Nt);
    B0=reshape(B,Nt,Nr);
    
 surfc((B0),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')    
axis tight
colorbar
%caxis([0 1])
%ylim([1 2600])
view([90 -90])

pause(0.3)
end