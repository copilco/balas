A=importdata('out1.txt');

B=reshape(A,1000,100);
figure
surf(B,'FaceColor','interp','EdgeColor','none')
 camlight left; lighting phong
 axis tight
 %axis([1 500 1 500 -1 1])
 %caxis([-0.7 0.7])
 
 
 
 colorbar('location','south')
 view(2)