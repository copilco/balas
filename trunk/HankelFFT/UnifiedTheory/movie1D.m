A=importdata('spout0.txt');

nr=600%120%1024
 %120

figure

    BB=reshape(A,nr,length(A)/nr);
 surfc((BB),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')    
axis tight
colorbar
caxis([0 1])
%ylim([1 2600])
view([90 90])

pause(0.3)
