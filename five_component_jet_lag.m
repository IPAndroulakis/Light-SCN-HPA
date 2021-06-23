%% Simulated individuals:shift-time for 5 components
load select_0215.txt
trans_days=zeros(length(select_0215(:,1)),5);
timeset=0:4800;
im = 12;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);
link=[0.25,1.0126 0.85];
for par=1:length(select_0215(:,1))
    link=select_0215(par,1:3);
A1=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
index=find(A1.x>2400);
t=A1.x(index);
indd=[1,8,9,16,19];
for yy=1:5
peak=[];
x1=A1.y(indd(yy),index);
x1=roundn(x1,-5);
    for bb=2:length(t)-1
        if x1(bb)>x1(bb-1) & x1(bb)>x1(bb+1)
            peak=[peak;t(bb)];
        end
    end
    record=[];
    phase0=mod(peak(1),24);
    mod(peak,24);
for jj=1:length(peak)
    if abs(mod(peak(jj),24)-phase0)>11.5
        record=[record;peak(jj)];
    end
end

if length(record)==0
   trans_days(yy)=0; 
else
trans_days(par,yy)=record(1)/24-100;
end
end
par
end
shifting_time_5comp=zeros(length(select_0215(:,1)),8);
shifting_time_5comp(:,1:3)=select_0215(:,1:3);
shifting_time_5comp(:,4:8)=trans_days;
writematrix(shifting_time_5comp)
%% plot
load shifting_time_5comp.txt
select=[];
for i=1:length(shifting_time_5comp(:,1))
        if shifting_time_5comp(i,4)>2
        if shifting_time_5comp(i,5)~=0
        if shifting_time_5comp(i,6)~=0
        if shifting_time_5comp(i,7)~=0
        if shifting_time_5comp(i,8)~=0
            select=[select;shifting_time_5comp(i,:)];
        end
        end
        end
        end
        end
end
figure
%X = categorical({'Per/Cry mRNA in core','VIP','Per/Cry mRNA in shell','AVP','CORT'});
%X = reordercats(X,{'Per/Cry mRNA in core','VIP','Per/Cry mRNA in shell','AVP','CORT'});

for i=1:5
    plot(i*ones(length(select(:,1)),1),select(:,i+3),'^','linewidth',2)
    hold on
end
m=mean(select(:,4:8))
s=std(select(:,4:8))
xi=[1 2 3 4 5];
plot(xi,m,'color',[0,0,0]+0.5,'linewidth',4,'Marker','o')
ylabel('Resynchronization time (days)')
names = {'Per (core)','VIP','Per (shell)','AVP','CORT'};
set(gca,'xtick',[1:5],'xticklabel',names)
set(gca,'fontsize', 14);
xlim([0 6])
ylim([0 105])


figure
for i=1:5
plot(m(i),s(i),'*','linewidth',2)
hold on
end
legend('Per/Cry (core)','VIP','Per/Cry (shell)','AVP','CORT')
ylabel('Std. dev.')
xlabel('Mean')
set(gca,'fontsize', 18);
%% ode function
function dzdt=gate_oscillator(t,z,pp,link)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

%permanant jet-lag
gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<2400)+((-1)*square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>=2400);

k=[0.381925841634067,3.06681115512865,0.349160115657042,4.38751198002070,0.456119328338128,1.00149445116453,0.848751479345601,0.803316775590302,0.724499581144526,0.180743285128964,2.99890752350971];
k2=[6.53671910362268,1.62794580266537,0.728122863970910];
h=[0.796500000000000;1.05770000000000;0.508400000000000;1.96270000000000;0.685700000000000;0.512900000000000;0.306900000000000;0.709700000000000;0.361800000000000;0.469500000000000;2.90000000000000;26.2000000000000;0.112403101000000;1.19876124000000;0.490000000000000;0.570000000000000;0.00329000000000000;0.0572000000000000;0.630000000000000;0.140100000000000;1.05770000000000;20;1.10000000000000;6;2;45;4;4;0.100000000000000;0.400000000000000;4;45;4;4;0.450000000000000;4;45;4;4;0.100000000000000;1.24000000000000;1.20000000000000;1.20000000000000;0.00300000000000000;15.7870000000000;81.6800000000000;0.140000000000000];
stress=1;

VIP=z(8);AVP=z(16);
CRH=z(17);ACTH=z(18);F_h=z(19);GRm_h=z(20);GR_h=z(21);FR_h=z(22);FRn_h=z(23);

v_l=link(1);v_c2=link(2);vcoe=link(3);


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
