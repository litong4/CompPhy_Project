%plot figures for Monte Carlo process in Ising model

clear; 

latt=20*20; 
num=latt*1000; 
x=[1.0,2.4]; 
h=0; 
for temp=x
    mag=zeros(2,num); 
    ene=zeros(2,num); 
    ll=0; 
    h=h+1; 
    for init=[0,1]
        ll=ll+1; 
        [mag(ll,:),ene(ll,:)]=textread(['..\benchmark\L20_init',num2str(init),'_temp',num2str(temp),'_mc.txt'],'%f %f\n','headerlines',1,'commentstyle','c++'); 
        figure(2*h-1); 
        plot((1:num)./latt,mag(ll,:)); 
        hold on; 
        figure(2*h); 
        plot((1:num)./latt,ene(ll,:)); 
        hold on; 
    end
    figure(2*h-1); 
    legend('Ordered initialization','Random initialization'); 
    xlabel('Number of Monte Carlo cycles'); 
    ylabel('Magnetization per spin'); 
    figure(2*h); 
    legend('Ordered initialization','Random initialization'); 
    xlabel('Number of Monte Carlo cycles'); 
    ylabel('Energy per spin'); 
end
