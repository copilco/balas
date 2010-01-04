%A=[1 0 0 0; 1  1 0 0; 1 1 1 0; 1 1 1 1]
%B=ones(4,3)
%C=A*B
%D=A'*B'

A=importdada('spout0.txt');
AA=A(:,1);

N=1024
MM=12
B=reshape(AA,N,MM);

surfc(B,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
axis tight
view([90 90]))

