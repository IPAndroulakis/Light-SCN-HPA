%% virtual population selection 
 t_extract=[2424 2428 2432 2436 2440 2444 2448];
 CORT_norminal=[0.2471    0.7116    3.1275    4.6603    3.5590    1.5939    0.2471];
 npop=10000;
p=sobolset(3);
R=[];
VarMin=[0, 0.8, 0.6];
VarMax=[1, 2 ,0.9];  
for i=1:npop
    r=p(i,:);
    r=VarMin+r.*(VarMax-VarMin);
    R=[R;r];
end
  
timeset=0:2400+1200;
im = 12;
B=1*ones(1,23);%Initial condition
select_0215=[];
option = odeset('RelTol', 1e-8);
for n=1:length(R(:,1))
    link=(R(n,:));
A=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
data=deval(A,t_extract);
CORT=data(19,:);
index=find(A.x>2400);
sse=sum(((CORT-CORT_norminal)).^2);

            cort=A.y(19,index);
            t_range=A.x(index);
            t=[];
            for i=2:length(index)-1
                if cort(i-1)<cort(i) & cort(i)>cort(i+1)
                    t(end+1)=t_range(i);
                end
            end
            delta=[];
            for j=1:length(t)-1
                delta(end+1)=t(j+1)-t(j); 
            end
period=mean(delta);
delta_mean=period-24;

 avp=A.y(16,index);
            t=[];
            for i=2:length(index)-1
                if avp(i-1)<avp(i) & avp(i)>avp(i+1)
                    t(end+1)=t_range(i);
                end
            end
            delta=[];
            for j=1:length(t)-1
                delta(end+1)=t(j+1)-t(j);
            end
period=mean(delta);
delta_mean2=period-24;

if sse < 2.8
   if abs(delta_mean2)<0.02 
   if abs(delta_mean)<0.02
select_0215=[select_0215;link]
   end
   end
end
end
writematrix(select_0215)
%% plot the selected population
load select_0215.txt
figure
%ttl={'Main view','Left view','Above view','Three-dimensional view'};
angle={[0,0],[-90,0],[0 90],[-37.5,30]};
for i=1:4
subplot(2,2,i)
    scatter3(select_0215(:,1),select_0215(:,2),select_0215(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.1)
   hold on
   scatter3(0.25,1.0126,0.85,200,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
    xlabel('v_l');ylabel('v_c_2');zlabel('v_c_o_e')
xlim([0 1]);ylim([0.5 2]);zlim([0.4 1]);
zticks([0.4 0.6 0.8 1])
set(gca,'fontsize', 18);
view(angle{i});%title(ttl{i});
end
%% plot CORT profiles in the simulated population
load select_0215.txt
select_0215_combined=select_0215;
nu=length(select_0215_combined(:,1))
pl=[];
figure
    scatter3(select_0215_combined(:,1),select_0215_combined(:,2),select_0215_combined(:,3))
    xlabel('vl');ylabel('vc2');zlabel('vcoe')
xlim([0 1]);ylim([0.5 2]);zlim([0.5 0.9]);
figure 
x2 = [0; 0;12;12];
y2 = [0 ;5;5;0];
P1=patch(x2,y2,'white','edgecolor','white');
x2 = [12;12;24;24];
y2 = [0 ;5;5;0];
P2=patch(x2,y2,[17 17 17]/255,'edgecolor','white');
alpha(0.3)
hold on
timeset=0:2400+1200;
im = 12;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);

for i=1:nu
    link=(select_0215_combined(i,1:3));
A=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
pl=[pl;A.x',A.y(19,:)'];
plot(A.x-2400,A.y(19,:),'color',[0,0,0]+0.5,'linewidth',2)
hold on
end
link=[0.25,1.0126 0.85];
A=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
plot(A.x-2400,A.y(19,:),'r-','linewidth',2) 
xlabel('Zeitgeber time (h)')
ylabel('CORT level')
xlim([0,24])
set(gca,'fontsize', 18);
%% Calculate resynchronization time
clc; clear all;
load select_0215.txt
timeset=0:4800;
im = 12;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);
shift_time=zeros(length(select_0215),1);
for co=1:length(select_0215(:,1))
A1=ode45(@(t,z)gate_oscillator(t,z,im,select_0215(co,:)),timeset,B,option);
index=find(A1.x>1200);
peak=[];
t=A1.x(index);
x1=A1.y(19,index);
    for bb=2:length(t)-1
        if x1(bb)>x1(bb-1) & x1(bb)>x1(bb+1)
            peak=[peak;t(bb)];
        end
    end
    phase0=mod(peak(1),24);
    record=[];
for jj=1:length(peak)
    if mod(peak(jj),24)<=(phase0-11.5)
        record=[record;peak(jj)];
    end
end
if length(record)==0
   shift_time(co)=0; 
else
shift_time(co)=record(1)/24-100;
end
shift_time(co)
end


select_0215_shift=zeros(length(select_0215(:,1)),4);
select_0215_shift(:,1)=select_0215(:,1);select_0215_shift(:,2)=select_0215(:,2);select_0215_shift(:,3)=select_0215(:,3);
select_0215_shift(:,4)=shift_time;
writematrix(select_0215_shift)
%% Alternating shift work
clc; clear all;
load select_0215.txt
load select_0215_shift.txt

timeset=0:2400+24*70;
im = 12;
B=1*ones(1,23);
option = odeset('RelTol', 1e-8);
average_gap_t=zeros(length(select_0215(:,1)),1);
for n=1:length(select_0215(:,1))
    link=(select_0215(n,1:3));
A1=ode45(@(t,z)gate_oscillator(t,z,im,link),timeset,B,option);
index=find(A1.x>2400);
per_vip=A1.y(19,index);
t_range=A1.x(index);
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
average_gap_t(n,1)=period-24;
end

select_0215_combined=zeros(length(select_0215(:,1)),5);
select_0215_combined(:,1:3)=select_0215(:,1:3);select_0215_combined(:,4)=select_0215_shift(:,4);select_0215_combined(:,5)=average_gap_t(:,1);
writematrix(select_0215_combined)
%% phase delay index test
clc; clear all
timeset=0:1200+24*200;
im = 12;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);

load select_0215_combined.txt
cee=zeros(length(select_0215_combined),1);
for len=1:length(select_0215_combined)
ee=0;
link=select_0215_combined(len,1:3);
for ls=1:1:12
A1=ode45(@(t,z)gate_oscillator_pdi(t,z,im,link,ls),timeset,B,option);
index=find(A1.x>1200);
per1=24;days1=199;
peak=zeros(1,days1-1);
t=A1.x(index)-1200;
x1=A1.y(1,index);
for td=1:days1-1
    t1=td*per1;
    t2=t1+(2*per1);
    k0=find(t>t1 & t<t2);
    teff=t(k0); 
    y=x1(1,k0);
    k=find(y<=max(y) & y>max(y)-(max(y)-mean(y))*0.1);
    ta=teff(k);
    peak(1,td)=ta(1,1)-td*per1;
end
for ii=1:length(peak)-1
    a(ii)=peak(1,ii+1)-peak(1,ii);
end
q=find(a>=10);
y=length(q);
if y==0
    ee=ee+1;
end
end
cee(len)=ee
end

select_0215_combinedd=zeros(length(select_0215_combined),6);
select_0215_combinedd(:,1:5)=select_0215_combined(:,1:5);
select_0215_combinedd(:,6)=cee;
writematrix(select_0215_combinedd)
 %% transient shift work
clc; clear all
timeset=0:2400+24*80;
im = 12;
B=1*ones(1,23);%Initial condition
option = odeset('RelTol', 1e-8);

load select_0215_combinedd.txt
del_phi=zeros(length(select_0215_combinedd(:,1)),1);
for oo=1:length(select_0215_combinedd(:,1))
    link=select_0215_combinedd(oo,1:3);
A0=ode45(@(t,z)gate_oscillator_ts(t,z,im,link,0),timeset,B,option);
A1=ode45(@(t,z)gate_oscillator_ts(t,z,im,link,1),timeset,B,option);
index0=find(A0.x>=2400);index1=find(A1.x>=2400);
cort0=A0.y(19,index0); cort1=A1.y(19,index1);
time0=A0.x(index0);time1=A1.x(index1);

peak0=[];
for nn=2:length(index0)-1
    if cort0(nn-1)<=cort0(nn)
        if cort0(nn)>=cort0(nn+1)
            peak0=[peak0;time0(nn)];
        end
    end
end
peak1=[];
for mm=2:length(index1)-1
    if cort1(mm-1)<=cort1(mm)
        if cort1(mm)>=cort1(mm+1)
            peak1=[peak1;time1(mm)];
        end
    end
end
del_phi(oo)=...
mod(max(abs(peak0(1:30)-peak1(1:30))),24)
end

select_0215_combineddd=zeros(length(select_0215_combinedd(:,1)),7);
select_0215_combineddd(:,1:6)=select_0215_combinedd;
select_0215_combineddd(:,7)=del_phi;
writematrix(select_0215_combineddd)
%% ode function
function dzdt=gate_oscillator(t,z,pp,link)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

%gate_output=square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;


%permanant jet-lag
gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<24*100)+((-1)*square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>=24*100);

%alternating shift work
%gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<24*50);
%if mod(t/24,7)<=5
%gate_output=square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;
%elseif mod(t/24,7)>5 
%   gate_output=(-1)*square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2;
%end

k=[0.381925841634067,3.06681115512865,0.349160115657042,4.38751198002070,0.456119328338128,1.00149445116453,0.848751479345601,0.803316775590302,0.724499581144526,0.180743285128964,2.99890752350971];
k2=[6.53671910362268,1.62794580266537,0.728122863970910]; % from rohit
h=[0.796500000000000;1.05770000000000;0.508400000000000;1.96270000000000;0.685700000000000;0.512900000000000;0.306900000000000;0.709700000000000;0.361800000000000;0.469500000000000;2.90000000000000;26.2000000000000;0.112403101000000;1.19876124000000;0.490000000000000;0.570000000000000;0.00329000000000000;0.0572000000000000;0.630000000000000;0.140100000000000;1.05770000000000;20;1.10000000000000;6;2;45;4;4;0.100000000000000;0.400000000000000;4;45;4;4;0.450000000000000;4;45;4;4;0.100000000000000;1.24000000000000;1.20000000000000;1.20000000000000;0.00300000000000000;15.7870000000000;81.6800000000000;0.140000000000000];
stress=1;

VIP=z(8);AVP=z(16);
CRH=z(17);ACTH=z(18);F_h=z(19);GRm_h=z(20);GR_h=z(21);FR_h=z(22);FRn_h=z(23);

k_avp=1;k_davp=1;
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

%% ode function for transient shift work
function dzdt=gate_oscillator_ts(t,z,pp,link,ind)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

day=8;
if ind==0
gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2);
else
gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<=2400)+((-1)*square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>2400 & t<2400+day*24)+(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>2400+day*24);
end    


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

    sc1*(v1be*(z(15)+v_c2*VIP^ce)/(k1be*(1+(z(11)/k1ie)^pe)+z(15)+v_c2*VIP^ce)-k1de*z(9));... 
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
%% ode function for phase delay index
function dzdt=gate_oscillator_pdi(t,z,pp,link,light_shift)

pin3=[0.0400    0.2500    2.0000    0.1000    1.0000 1.5 1.5 1.0126    3.0919 1.1 1.1 1.05 1.1 1.1 1.1 1.1 1.02 1.036]; %original delta=3.6h, entraiment range=20h
lis=pin3(1);v_l=pin3(2);K_l=pin3(3); v_c1=pin3(4); cm=pin3(5); k_vip=pin3(6);k_dvip=pin3(7);
v_c2=pin3(8); ce=pin3(9); sc1=pin3(10); sc2=pin3(11);sc3=pin3(12); sc4=pin3(13);sc5=pin3(14);sc6=pin3(15); sc7=pin3(16);
f=pin3(17);g=pin3(18);

v1bo=9; k1bo=1; k1io=0.56;k1do=0.12;k2bo=0.3;k2do=0.05;k3do=0.12;k2to=0.24;k3to=0.02;v4bo=3.6;k4bo=2.16;k4do=0.75;k5bo=0.24;k5do=0.06;k6do=0.12;k5to=0.45;k6to=0.06;k6po=0.09;k7po=0.003;k7do=0.009;po=8;qo=2;ro=3; %o refers to original value
v1bm=v1bo*f; k1bm=k1bo*f; k1im=k1io*f;k1dm=k1do*f;k2bm=k2bo*f;k2dm=k2do*f;k3dm=k3do*f;k2tm=k2to*f;k3tm=k3to*f;v4bm=v4bo*f;k4bm=k4bo*f;k4dm=k4do*f;k5bm=k5bo*f;k5dm=k5do*f;k6dm=k6do*f;k5tm=k5to*f;k6tm=k6to*f;k6pm=k6po*f;k7pm=k7po*f;k7dm=k7do*f;pm=po;qm=qo;rm=ro; %f refers to scaled factor
v1be=v1bo*g; k1be=k1bo*g; k1ie=k1io*g;k1de=k1do*g;k2be=k2bo*g;k2de=k2do*g;k3de=k3do*g;k2te=k2to*g;k3te=k3to*g;v4be=v4bo*g;k4be=k4bo*g;k4de=k4do*g;k5be=k5bo*g;k5de=k5do*g;k6de=k6do*g;k5te=k5to*g;k6te=k6to*g;k6pe=k6po*g;k7pe=k7po*g;k7de=k7do*g;pe=po;qe=qo;re=ro; %g refers to scaled factor

gate_output=(square(t*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t<1200)+(square((t-light_shift)*(2*pi)/24,(pp/24)*100)*lis/2+lis/2).*(t>1200);


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