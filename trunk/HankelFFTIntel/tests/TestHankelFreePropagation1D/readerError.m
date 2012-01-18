A=importdata('normrecord.txt');
figure
plot(A,'LineWidth',4)
title('Error Norm R=200, Nr=1000, dt=0.01','FontSize',13)
xlabel('Iterations','FontSize',13)
set(gca,'FontSize',13)
