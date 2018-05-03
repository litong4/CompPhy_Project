%plot figures for the phase transition in Ising model

clear; 

for size=[20 40 60 80 100 140]
    v=zeros(7,6); 
    ll=0; 
    x=2.1:0.02:2.4;
    for temp=x
        ll=ll+1; 
        v(ll,:)=transpose(textread(['..\benchmark\tran_size',num2str(size),'_temp',num2str(temp),'_sum.txt'],'%*s %f\n','delimiter',':','headerlines',1)); 
        v(ll,:)=v(ll,:)./(size*size); 
    end

    figure(1); 
    plot(x,v(:,1),'.-'); %magnetization
    hold on; 
    figure(2); 
    plot(x,v(:,3),'.-'); %energy
    hold on; 
    figure(3); 
    plot(x,v(:,5),'.-'); %C_v
    hold on; 
    figure(4); 
    plot(x,v(:,6),'.-'); %susceptiblity
    hold on; 
end 

figure(1); 
legend('L=20','L=40','L=60','L=80','L=100','L=140','Location','northeast'); 
xlabel('Temperature'); 
ylabel('Magnetization per spin'); 
figure(2); 
legend('L=20','L=40','L=60','L=80','L=100','L=140','Location','northwest'); 
xlabel('Temperature'); 
ylabel('Energy per spin'); 
figure(3); 
legend('L=20','L=40','L=60','L=80','L=100','L=140','Location','northwest');
xlabel('Temperature'); 
ylabel('Heat Capacity per spin'); 
figure(4); 
legend('L=20','L=40','L=60','L=80','L=100','L=140','Location','northwest');
xlabel('Temperature'); 
ylabel('Susceptibility per spin'); 
