a=importdata('out0.txt');

N=50
aa=reshape(a,N,N);

figure
surfc(aa','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
axis tight
%caxis([0 1])
%ylim([1 2600])
view([90 90])



b=importdata('out1.txt');

N=50
bb=reshape(b,N,N);

figure
surfc(bb','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
axis tight
%caxis([0 1])
%ylim([1 2600])
view([90 90])


