%%% Simple plot for energy

ene=importdata('outEne.txt');

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
    'Color','w');
plot(ene,'o');
xlabel('Temporal iterations','fontsize',12)
ylabel('Energy (a.u.)','fontsize',12)
title('Energy of ground state','fontsize',16)
grid on