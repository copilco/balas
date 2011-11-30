A=importdata('out2.txt');

figure
N=1000
for i=1:100
    A1=A((i-1)*N+1:i*N,:);
    
    hold off
    plot(A1(:,2),A1(:,3),'ro')
    hold on
    plot(A1(:,2),A1(:,3),'g*')
    
    pause(0.05)
end