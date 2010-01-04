figure
a0=importdata('fileinX.txt');
a1=importdata('fileinQ.txt');
a2=importdata('fileoutX.txt');

b0=reshape(a0,60,40);
b1=reshape(a1,60,40);
%bb1=fftshift(b1);
b2=reshape(a2,60,40);

const=1e-3;
const2=1e-3;

subplot (1,3,1)
p=surfc(b0)
set(p,'FaceColor','red','EdgeColor','none');
%daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
title('Original')

subplot (1,3,2)
p1 = surfc(b1);
set(p1,'FaceColor','green','EdgeColor','none');
%daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
title('After FFT in 1d Forward')

subplot (1,3,3)
p2 = surfc(b2);
set(p2,'FaceColor','red','EdgeColor','none');
%daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

title('After FFT in 1d Back')