rmax=xlsread('proj2_data.xlsx','one','A2:A6');
energy=xlsread('proj2_data.xlsx','one','B2:D6');
figure(1);
plot(rmax,energy(:,1),'b.-'); 
hold on; 
plot(rmax,energy(:,2),'r.-'); 
plot(rmax,energy(:,3),'k.-'); 
plot([1 5],[3.0 3.0],'b--'); 
plot([1 5],[7.0 7.0],'r--'); 
plot([1 5],[11.0 11.0],'k--');
axis([-inf,inf,0,40]); 
xlabel('$\rho_{\max}$','FontName', 'Times New Roman','Interpreter', 'LaTeX');
ylabel('$\lambda$','FontName', 'Times New Roman','Interpreter', 'LaTeX');
legend('Ground state','1st excited state','2nd excited state');