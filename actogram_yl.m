function actogram_yl(t,x1,th1,per1,days1,rep1,x2,th2,per2,days2,rep2,x3,th3)
    
%%% Plot the actogram of the coupled oscillator model (coupled_Model.m) 
% t: Time
% x1:Time series 1 for M activity
% th1: Threshold for defining M activity
% per1 = period (typically 24h)
% days1: Number of days for M activity
% rep1: Number of repetition for M activity (1 for single actogram, 2 for double actogram)
% x2:Time series 2 for E activity
% th2: Threshold for defining E activity
% per2 = period (typically 24h)
% days2: Number of days for E activity
% rep2: Number of repetition for E activity (1 for single actogram, 2 for double actogram)

clf;
peak=zeros(3,days1-1);
%%% Plot data
figure
for td=1:days1-1
    t1=td*per1;
    t2=t1+(rep1*per1);
    k0=find(t>t1 & t<t2);
    teff=t(k0); 
    y=x1(k0);
    k=find(y<=max(y) & y>max(y)-(max(y)-mean(y))*0.1);
    ta=teff(k);
    peak(1,td)=ta(1,1)-td*per1;

    h=days1-td;

    scatter(ta-td*per1,h*ones(1,length(k)),8,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.1)
    hold on;
    
end
%**********************************************************************

for td=1:days2-1     
    t1=td*per2;
    t2=t1+(rep2*per2);
    k=find(t>t1 & t<t2);
    teff=t(k); 
    y=x2(k);
    k=find(y<=max(y) & y>max(y)-(max(y)-mean(y))*0.1);
   ta=teff(k);
    peak(2,td)=ta(1,1)-td*per2;

    h=days2-td;

     scatter(ta-td*per2,h*ones(1,length(k)),8,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.1)
    hold on;
    
end

days3=days2; per3=per2; rep3=rep2;
for td=1:days3-1
        
    t1=td*per3;
    t2=t1+(rep3*per3);
    k=find(t>t1 & t<t2);
    teff=t(k);
    y=x3(k);
    k=find(y<=max(y) & y>max(y)-(max(y)-mean(y))*0.1);
   ta=teff(k);
 peak(3,td)=ta(1,1)-td*per3;
 
    h=days3-td;

    %p3=plot(ta-td*per3,h*ones(1,length(k)),'k.','Markersize',8);
    %p3.Color(4) = 0.5;
    scatter(ta-td*per3,h*ones(1,length(k)),8,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.1); 
    hold on;
    
end

%*********************************************************************

%%% Axis limits and labels

xlim([0 48])
ylim([0 days1])
xticks([0 6 12 18 24 30 36 42 48])
xticklabels({'0','6','12','18','0','6','12','18','0'})
set(gca,'YTickLabel',{''})

xlabel('Peaking time (h)','fontsize',18);
ylabel('Days','fontsize',18);

ylabh = get(gca,'yLabel'); 
set(ylabh,'Position',get(ylabh,'Position') - [2.5 0 0]);  % shift y-label to left 

figure
plot((1:days1-1),peak(1,:),'bo','linewidth',2)
hold on
plot((1:days1-1),peak(2,:),'ro','linewidth',2)
plot((1:days1-1),peak(3,:),'ko','linewidth',2)
ylabel('Peaking time (h)','fontsize',18);
xlabel('Days','fontsize',18);
