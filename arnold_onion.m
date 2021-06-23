%% arnold onion
timeset=0:2400+24*7;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);
% link=[0.25,1.0126,0.85];
link=[0.25,0,0.85;0.25,0.5,0.85;0.25,1,0.85;0.25,1.5,0.85;0.25,2,0.85;...
0.25,1.0126,0;0.25,1.0126,0.4;0.25,1.0126,0.8;0.25,1.0126,1.2;0.25,1.0126,1.6];
deltaaa=[];
for hh=1:length(link(:,1))
    hh
    link_p=link(hh,:);
geo=[];geo1=[];geo2=[];
%for E_T=18:0.01:27 %for drawing onions
%    for pp=0:0.02:1 %for drawing onions
for E_T=18:0.2:27 % for geting parameter functions 
    for pp=0:0.1:1 % for geting parameter functions 
        E(1)=pp;E(2)=E_T;
        A=ode45(@(t,z)gate_oscillator(t,z,link_p,E),timeset,B,option);
         index=find(A.x>=2400);
per_vip=A.y(19,index);
t_range=A.x(index);
t=[];
for i=2:length(index)-1
    if per_vip(i-1)< per_vip(i) && per_vip(i)> per_vip(i+1)
        t(end+1)=t_range(i);
    end
end
delta=[];
for j=1:length(t)-1
    delta(end+1)=t(j+1)-t(j);
end
period=mean(delta);
        if abs(period-E_T)<0.1
            geo=[geo;E_T, pp];
        end
        %
per_vip=A.y(8,index);
t_range=A.x(index);
t=[];
for i=2:length(index)-1
    if per_vip(i-1)< per_vip(i) && per_vip(i)> per_vip(i+1)
        t(end+1)=t_range(i);
    end
end
delta=[];
for j=1:length(t)-1
    delta(end+1)=t(j+1)-t(j);
end
period=mean(delta);
        if abs(period-E_T)<0.1
            geo1=[geo1;E_T, pp]
        end
        %
        per_vip=A.y(16,index);
t_range=A.x(index);
t=[];
for i=2:length(index)-1
    if per_vip(i-1)< per_vip(i) && per_vip(i)> per_vip(i+1)
        t(end+1)=t_range(i);
    end
end
delta=[];
for j=1:length(t)-1
    delta(end+1)=t(j+1)-t(j);
end
period=mean(delta);
        if abs(period-E_T)<0.1
            geo2=[geo2;E_T, pp];
        end
        
 index=[];            
    end
end
deltaaa=[deltaaa;length(geo1(:,1))-length(geo2(:,1)),length(geo2(:,1))-length(geo(:,1))]
end
%% 
figure
subplot(1,2,1)
x=link(1:5,2)
y=deltaaa(1:5,1)
plot(x,y,'ks','linewidth',2)
hold on
plot(x,y,'k-','linewidth',2)
xlabel('v_c_2')
ylabel('Area difference between VIP and AVP')
set(gca,'fontsize', 16);

subplot(1,2,2)
x=link(6:10,3)
y=deltaaa(6:10,2)
plot(x,y,'ks','linewidth',2)
hold on
plot(x,y,'k-','linewidth',2)
xlabel('v_c_o_e')
ylabel('Area difference between AVP and CORT')
set(gca,'fontsize', 18);
%%
geoo1=[];geoo=[];
for i=1:length(geo1(:,1))
    if geo1(i,1)>18.4
        geoo1=[geoo1;geo1(i,:)];
    end
end
for i=1:length(geo(:,1))
    if geo(i,1)>20.66
        geoo=[geoo;geo(i,:)];
    end
end
figure
set(gca,'fontsize', 18);
subplot(2,2,1)
plot(geoo1(:,1),geoo1(:,2),'b.','linewidth',1)
set(gca,'fontsize', 18);
xlim([18,28])
ylim([0 1]);title('VIP')
xlabel('T (h)')
ylabel('Photoperiod x')
subplot(2,2,2)
plot(geo2(:,1),geo2(:,2),'r.','linewidth',1)
set(gca,'fontsize', 18);
xlim([18,28])
ylim([0 1]);title('AVP')
xlabel('T (h)')
ylabel('Photoperiod x')
subplot(2,2,3)
plot(geoo(:,1),geoo(:,2),'k.','linewidth',1)
set(gca,'fontsize', 18);
xlim([18,28])
ylim([0 1]);title('CORT')
xlabel('T (h)')
ylabel('Photoperiod x')
subplot(2,2,4)
plot(geoo1(:,1),geoo1(:,2),'b.','linewidth',1)
hold on
plot(geo2(:,1),geo2(:,2),'r.','linewidth',1)
alpha(0.5);
plot(geoo(:,1),geoo(:,2),'k.','linewidth',1)
alpha(0.5);
xlim([18,28])
ylim([0 1]);title('combine')
xlabel('T (h)')
ylabel('Photoperiod x')
set(gca,'fontsize', 18);

%% ode function
function dzdt=gate_oscillator(t,z,link,E)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

%lis=0.1;
pp=E(1);T=E(2);
gate_output=square(t*(2*pi)/T,pp*100)*lis/2+lis/2;

k=[0.381925841634067,3.06681115512865,0.349160115657042,4.38751198002070,0.456119328338128,1.00149445116453,0.848751479345601,0.803316775590302,0.724499581144526,0.180743285128964,2.99890752350971];
k2=[6.53671910362268,1.62794580266537,0.728122863970910]; % from rohit
h=[0.796500000000000;1.05770000000000;0.508400000000000;1.96270000000000;0.685700000000000;0.512900000000000;0.306900000000000;0.709700000000000;0.361800000000000;0.469500000000000;2.90000000000000;26.2000000000000;0.112403101000000;1.19876124000000;0.490000000000000;0.570000000000000;0.00329000000000000;0.0572000000000000;0.630000000000000;0.140100000000000;1.05770000000000;20;1.10000000000000;6;2;45;4;4;0.100000000000000;0.400000000000000;4;45;4;4;0.450000000000000;4;45;4;4;0.100000000000000;1.24000000000000;1.20000000000000;1.20000000000000;0.00300000000000000;15.7870000000000;81.6800000000000;0.140000000000000];
stress=1;

VIP=z(8);AVP=z(16);
CRH=z(17);ACTH=z(18);F_h=z(19);GRm_h=z(20);GR_h=z(21);FR_h=z(22);FRn_h=z(23);
%v_l=0;v_l=0.25;
%v_c2=0;v_c2=1.0126;
%vcoe=0;vcoe=0.85;

k_avp=1;k_davp=1;
v_l=link(1);v_c2=link(2);vcoe=link(3);
%k_vip=link(1);k_avp=link(2);vcoe=link(3);

%v_l=0.*(t<1200+24*300)+(0.25).*(t>1200+24*300);
%v_c2=0.*(t<1200+24*100)+(1.0126).*(t>1200+24*100);
%vcoe=0.*(t<1200+24*200)+(0.85).*(t>1200+24*200);

%v_l=0.25.*(t<1200+24*100)+(0).*(t>1200+24*100);
%v_c2=1.0126.*(t<1200+24*300)+(0).*(t>1200+24*300);
%vcoe=0.85.*(t<1200+24*200)+(0).*(t>1200+24*200);

dzdt=[v1bm*(z(7)+v_c1*VIP^cm)/(k1bm*(1+(z(3)/k1im)^pm)+z(7)+v_c1*VIP^cm)-k1dm*z(1)+gate_output;...%lis,v_l,K_l,v_c1,K_c1 
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
        k_avp*z(10)-k_davp*AVP;...

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