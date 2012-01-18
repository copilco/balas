A=importdata('out1.txt');

figure
N=1000
for i=1:100
    A1=A((i-1)*N+1:i*N,:);
    
    hold off
    plot(A1(:,1),A1(:,2),'ro')
    hold on
    plot(A1(:,3),A1(:,4),'g*')
    axis tight
    ylim([0 0.03])
    
    pause(0.05)
end