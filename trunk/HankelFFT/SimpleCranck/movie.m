A=importdata('movie.txt');

B=A;
BB=reshape(B,200,length(B)/200);

figure
surf(BB,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

%'FaceColor','interp',...
%          'EdgeColor','none')
 grid off
 view([0 90])