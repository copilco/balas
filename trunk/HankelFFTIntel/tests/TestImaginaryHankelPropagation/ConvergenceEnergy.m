%%% Simple plot for energy

ene = importdata('outEne.txt');
proj = importdata('outProj.txt');

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
    'Color','w');
plot(ene(:,1),'ro');
xlabel('Temporal iterations','fontsize',12)
ylabel('Energy (a.u.)','fontsize',12)
title('Energy of ground state','fontsize',16)
grid on
hleg1 = legend('Energy','Logaritmic error in energy');

%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
    'Color','w');
plot(ene(:,2),'x');
xlabel('Temporal iterations','fontsize',12)
ylabel('Energy (a.u.)','fontsize',12)
title('Logaritmic error of the energy of ground state','fontsize',16)
grid on



%%

MAX = max(proj);
MIN = min(proj);

Band = MAX-MIN;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.7],...
    'Color','w');
plot(proj(:,1),'go')
xlabel('Temporal iterations','fontsize',12)
ylabel('Error (a.u.)','fontsize',12)
title('Error of the projection between the two ground states','fontsize',16)
grid on

