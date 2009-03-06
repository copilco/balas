A=importdata('pout4.txt');
n=100
for i=1:40
    B=A((i-1)*n*n+1:i*n*n);
    C=reshape(B,n,n);
    %figure
surf(C,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

%daspect([5 5 1])
axis tight
view(2)
pause(0.5)
end