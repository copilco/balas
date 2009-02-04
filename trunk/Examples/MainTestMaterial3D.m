A=importdata('out4.txt');

n=20
for i=1:25
    clf
    B=A((i-1)*n^3+1:i*n^3);
    C=reshape(B,n,n,n);
    
%    [x,y,z,v] = flow;
subplot(3,3,[1 2 4 5])
p = patch(isosurface(C,0.5e-4));
isonormals(C,p)
set(p,'FaceColor','red','EdgeColor','none');
%daspect([1 1 1])
xlim([1 n])
ylim([1 n])
zlim([1 n])
view(3);% axis tight
camlight 
lighting gouraud
alpha(0.7)


 d1=sum(C,1);
 f1(:,:)=d1(1,:,:);
 subplot(3,3,[3 6])
 surfc(f1, 'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
%daspect([5 5 1])
axis tight
view(2)
%camlight left

d2=sum(C,2);

 f2(:,:)=d2(:,1,:);
g2=sum(f2,2);
%h2(:)=g2(:,1);

 
 subplot(3,3,[7 8])
 plot(g2(:,1))
 axis tight
% surfc(f2, 'FaceColor','interp',...
% 'EdgeColor','none',...
% 'FaceLighting','phong')
%daspect([5 5 1])
%axis tight
%view(2)
%camlight left

pause(0.8)
end