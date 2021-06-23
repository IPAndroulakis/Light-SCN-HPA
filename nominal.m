%% Dynamic profiles of clock gene mRNAs, VIP, AVP and CORT
clear all;clc;close all
link=[0.25,1.0126,0.85]; %three coupling strengths in nominal parameter set
timeset=0:2400+24*20;
B=1*ones(1,23); 
option = odeset('RelTol', 1e-8);
im = 12; %photoperiod
A=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);

in_ss=find(A.x>=1200);
time_ss=A.x(in_ss);
CORT_ss=A.y(19,in_ss);
delta=max(CORT_ss)-min(CORT_ss); %amplitude of CORT

figure
time=A.x;VIP=A.y(8,:);AVP=A.y(16,:);Per_vip=A.y(1,:);Per_avp=A.y(9,:);Bmal1_vip=A.y(4,:);Bmal1_avp=A.y(12,:);CORT=A.y(19,:);
x2 = [0; 0;12;12;24; 24;36;36];
y2 = [0 ;7;7;0;0 ;7;7;0];
P1=patch(x2,y2,'white','edgecolor','white');
x2 = [12;12;24;24;36;36;48;48];
y2 = [0 ;7;7;0;0 ;7;7;0];
P2=patch(x2,y2,[17 17 17]/255,'edgecolor','white');
alpha(0.3)
hold on

linew=4;
plot(time-2400,VIP,'y','LineWidth',linew)
hold on
plot(time-2400,AVP,'g','LineWidth',linew)
plot(time-2400,Per_vip,'r-','LineWidth',linew)
plot(time-2400,Per_avp,'r--','LineWidth',linew)
plot(time-2400,Bmal1_vip,'b-','LineWidth',linew)
plot(time-2400,Bmal1_avp,'b--','LineWidth',linew)
plot(time-2400,CORT,'k','LineWidth',linew)
legend('Day','Night','VIP','AVP','Per mRNA (core)','Per mRNA (shell)','Bmal1 mRNA (core)','Bmal1 mRNA (shell)','CORT')
ylabel('Concentrations')
xlabel('Zeitgeber time(h)')
xlim([0,48])
ylim([0,6])
set(gca,'fontsize', 18);
%% Synchronization dynamics & jat leg
timeset=0:1200+24*500;
im = 12;
B=1*ones(1,23);
option = odeset('RelTol', 1e-8);
link=[0.25,1.0126 0.85];
A1=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
figure
plot(A1.x-2400,A1.y(1,:),'linewidth',2)
hold on
plot(A1.x-2400,A1.y(9,:),'linewidth',2)
plot(A1.x-2400,A1.y(19,:),'linewidth',2)
xlim([0,24*20])

index=find(A1.x>1200);
time=A1.x(index);
x1=A1.y(8,index); %VIP
th1=1.7*mean(x1);
per1=24;days1=400;rep1=2;
x2=A1.y(16,index); %AVP
th2=1.5*mean(x2);
per2=24;days2=400;rep2=2;
x3=A1.y(19,index); %CORT
th3=1.21*mean(x3');
actogram_yl(time-1200,x1,th1,per1,days1,rep1,x2,th2,per2,days2,rep2,x3,th3)
%% shift work
timeset=0:24*7*20+24*7*80;
im = 12;
B=1*ones(1,23);
option = odeset('RelTol', 1e-8);
link=[0.25,1.0126 0.85]; % uncomment for simulating group A
%link=[0.6875,1.775,0.84375]; %uncomment for simulating group B
A1=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
index=find(A1.x>24*7*20);
time=A1.x(index);
x1=A1.y(8,index);
th1=1.7*mean(x1);
per1=24;days1=400;rep1=2;
x2=A1.y(16,index);
th2=1.5*mean(x2);
per2=24;days2=400;rep2=2;
x3=A1.y(19,index);
th3=1.21*mean(x3');
actogram_yl(time-24*7*20,x1,th1,per1,days1,rep1,x2,th2,per2,days2,rep2,x3,th3)

%% ode function
function dzdt=gate_oscillator(t,z,pp,link)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

%normal L12/D12
gate_output=square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2; 

%jet-lag with 6 h phase advance 
%gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<1200+24*200)+(square((t+6)*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>1200+24*200);

%shift work
%gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<24*7*20+24*7*10);
%if t>= 24*7*20+24*7*10 & mod(t/24,7)<=5
%gate_output=square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;
%elseif t>= 24*7*20+24*7*10 & mod(t/24,7)>5 
%   gate_output=(-1)*square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;
%end

% k, k2, h parameter values are taken from (Rao & Androulakis, 2019) 
k=[0.381925841634067,3.06681115512865,0.349160115657042,4.38751198002070,0.456119328338128,1.00149445116453,0.848751479345601,0.803316775590302,0.724499581144526,0.180743285128964,2.99890752350971];
k2=[6.53671910362268,1.62794580266537,0.728122863970910]; 
h=[0.796500000000000;1.05770000000000;0.508400000000000;1.96270000000000;0.685700000000000;0.512900000000000;0.306900000000000;0.709700000000000;0.361800000000000;0.469500000000000;2.90000000000000;26.2000000000000;0.112403101000000;1.19876124000000;0.490000000000000;0.570000000000000;0.00329000000000000;0.0572000000000000;0.630000000000000;0.140100000000000;1.05770000000000;20;1.10000000000000;6;2;45;4;4;0.100000000000000;0.400000000000000;4;45;4;4;0.450000000000000;4;45;4;4;0.100000000000000;1.24000000000000;1.20000000000000;1.20000000000000;0.00300000000000000;15.7870000000000;81.6800000000000;0.140000000000000];
stress=1;

VIP=z(8);AVP=z(16);
CRH=z(17);ACTH=z(18);F_h=z(19);GRm_h=z(20);GR_h=z(21);FR_h=z(22);FRn_h=z(23);

v_l=link(1);v_c2=link(2);vcoe=link(3);

%uncomment for simulating synchronization dynamics
%v_l=0.*(t<1200+24*300)+(0.25).*(t>1200+24*300);
%v_c2=0.*(t<1200+24*100)+(1.0126).*(t>1200+24*100);
%vcoe=0.*(t<1200+24*200)+(0.85).*(t>1200+24*200);

dzdt=[v1bm*(z(7)+v_c1*VIP^cm)/(k1bm*(1+(z(3)/k1im)^pm)+z(7)+v_c1*VIP^cm)-k1dm*z(1)+v_l*gate_output/(gate_output+K_l);...%lis,v_l,K_l,v_c1,K_c1 
    k2bm*z(1)^qm-k2dm*z(2)-k2tm*z(2)+k3tm*z(3);...
    k2tm*z(2)-k3tm*z(3)-k3dm*z(3);...
    v4bm*z(3)^rm/(k4bm+z(3)^rm)-k4dm*z(4);...
    k5bm*z(4)-k5dm*z(5)-k5tm*z(5)+k6tm*z(6);...
    k5tm*z(5)-k6tm*z(6)-k6dm*z(6)+k7pm*z(7)-k6pm*z(6);...
    k6pm*z(6)-k7pm*z(7)-k7dm*z(7);...
        k_vip*z(2)-k_dvip*VIP;...

    sc1*(v1be*(z(15)+v_c2*VIP^ce)/(k1be*(1+(z(11)/k1ie)^pe)+z(15)+v_c2*VIP^ce)-k1de*z(9));... % v1be~9, k1be~1, k1ie~0.56, pe~8, k1de~0.12
    sc2*(k2be*z(9)^qe-k2de*z(10)-k2te*z(10)+k3te*z(11));...
    sc3*(k2te*z(10)-k3te*z(11)-k3de*z(11));...
    sc4*(v4be*z(11)^re/(k4be+z(11)^re)-k4de*z(12));...
    sc5*(k5be*z(12)-k5de*z(13)-k5te*z(13)+k6te*z(14));...
    sc6*(k5te*z(13)-k6te*z(14)-k6de*z(14)+k7pe*z(15)-k6pe*z(14));...
    sc7*(k6pe*z(14)-k7pe*z(15)-k7de*z(15));...
        z(10)-AVP;...

% Corticotropin Releasing Hormone
k(1)*stress*k2(1)/(k2(1)+FRn_h)-k(3)*CRH*(1+vcoe*AVP.^3/(1+AVP.^3))/(k(4)+CRH);...
% ACTH
k(5)*k2(2)*CRH/(k2(2)+FRn_h^1)-1*k(6)*ACTH/(1*k(7)+ACTH);... 
% Corticosterone
1*k2(3)*ACTH-k(9)*F_h/(1*k(10)+F_h);...
% GR mNRA in HPA
h(11)*(1.0)*(1-FRn_h/(h(12)+FRn_h))-h(13)*GRm_h;... 
% GR Protein in HPA
h(14)*GRm_h+h(15)*h(16)*FRn_h-h(17)*(F_h)*GR_h-h(18)*GR_h;... 
% Cytoplasmic CST-Bound Receptor in HPA
h(17)*(F_h)*GR_h-h(19)*FR_h;... 
% Nuclear CST-Bound Receptor
h(19)*(FR_h)-h(16)*FRn_h];
end
