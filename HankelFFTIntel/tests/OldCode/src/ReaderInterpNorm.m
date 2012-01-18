A=importdata('out2.txt');


%%
plot(A(:,1),log10(abs(A(:,2))),'r*')
xlabel('Number of points Nr')
ylabel('Error (log10 scale)')
title('Error in normalization of Hankel Interpolated')