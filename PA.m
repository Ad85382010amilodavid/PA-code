A=5;%recriutment rate 
mu=0.03;%natural death coefficient
d=0.02;
beta_1=0.8;%transmission infectivity rate
beta_2=0.9;%death rate  
r_1=0.1;% gama recovery induced rate
r_2=0.1;% rate at which recovered return to suceptibles
e_1=0.2;% delta The rate at which addictive develop alcoholic cardiomyopathy disease
e_2=0.2;% Rate at which alcoholic cardiomyopathy diseased receive treatment
c_1=0.6;
c_2=0.6;
s_1=0.4;
s_2=0.4;
z=0.03;
k_1=0.6;
k_2=0.4;
gama=0.002;
 

%fdefun=@(t,y)[a-b*y(2)-u*y(1); -q*y(2)-p*y(2)*y(3); q*y(2)*y(3)-(u+P+d)*y(3); p*y(3)-(u+g)*y(4); (g*y(4))-(d+g)*y(5)];
fdefun=@(t,y)[A+z*y(2)+(r_1)*y(7)-(beta_2)*y(3)-mu*y(1)-gama*y(1); gama*y(1)-z*y(2)-mu*y(2)-(beta_1)*y(4)+(r_2)*y(7); (beta_2)*y(3)-(mu+d)*y(3)-(c_1)*y(3)-(s_1)*y(3);(beta_1)*y(4)-(mu+d)*y(4)-(c_2)*y(4)-(s_2)*y(4); (c_1)*y(3)+(k_2)*y(6)-(k_1)*y(5)-(mu+d)*y(5)+(c_2)*y(4); (s_1)*y(3)+(k_1)*y(5)-(k_2)*y(6)-(mu+d)*y(6)+(s_2)*y(4);(e_1)*y(5)+(e_2)*y(6)-(r_1)*y(7)-(r_2)*y(7)-mu*y(7)];

alpha1=0.4; 
alpha=0.6;
alpha2=1;
t0 = 0 ; tfinal =400 ; y0 = [3828 ; 1914; 168; 336;672;338;536] ;
h = 2^(-6) ;
[t1, y_fde121] = fde12(alpha1,fdefun,t0,tfinal,y0, h);
[t, y_fde12] = fde12(alpha,fdefun,t0,tfinal,y0, h);
[t2, y_fde122] = fde12(alpha2,fdefun,t0,tfinal,y0, h);

 
figure(1)
plot(t,y_fde121(1,:),'o-b',t, y_fde12(1,:),'b',t, y_fde122(1,:),'--b',t,y_fde121(2,:),'o-m',t, y_fde12(2,:),'m',t, y_fde122(2,:),'--m',t,y_fde121(3,:),'o-y',t, y_fde12(3,:),'y',t, y_fde122(3,:),'--y',t,y_fde121(4,:),'o-r',t, y_fde12(4,:),'r',t, y_fde122(4,:),'--r',t,y_fde121(5,:),'o-g',t, y_fde12(5,:),'g',t, y_fde122(5,:),'--g',t,y_fde121(6,:),'o-g',t, y_fde12(6,:),'g',t, y_fde122(6,:),'-g',t,y_fde121(7,:),'o-c',t, y_fde12(7,:),'c',t, y_fde122(7,:),'--c');
legend('S_C(t)alpha=0.4','S_C(t)alpha=0.6','S_C(t)alpha=1','S_H(t)alpha=0.4','S_H(t)alpha=0.6','S_H(t)alpha=1','I_C(t)alpha=0.4','I_C(t)alpha=0.6','I_C(t)alpha=1','I_H(t)alpha=0.4','I_H(t)alpha=0.6','I_H(t)alpha=1','T_A_1(t)alpha=0.4','T_A_1(t)alpha=0.6','T_A_1(t)alpha=1','T_A_2(t)alpha=0.4','T_A_2(t)alpha=0.6','T_A_2(t)alpha=1','R(t)alpha=0.4','R(t)alpha=0.6','R(t)alpha=1');
xlabel('time')
ylabel('population')

%%Saving
%set(gca, 'YScale', 'log')
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'ZeroToOne', 'pdf')



