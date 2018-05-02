%plot figures for Monte Carlo process in Ising model

clear; 

num=20*20*1000; 
x=[1.0,2.4]; 
h=0; 
for temp=x
    mag=zeros(2,num); 
    ene=zeros(2,num); 
    ll=0; 
    h=h+1; 
    for init=[1]
        ll=ll+1; 
        [mag(ll,:),ene(ll,:)]=textread(['..\benchmark\L20_init',num2str(init),'_temp',num2str(temp),'_mc.txt'],'%f %f\n','headerlines',1,'commentstyle','c++'); 
        figure(2*h-1); 
        histogram(mag(ll,:),100); 
        figure(2*h); 
        histogram(ene(ll,:),100); 
    end
end
