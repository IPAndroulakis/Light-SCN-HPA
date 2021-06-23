%% Plot resynchronization time in parametric space
load select_0215_shift.txt
select_0215_shift1=[];
for n=1:length(select_0215_shift(:,1))
    if select_0215_shift(n,4)<100  & select_0215_shift(n,4)>0
  select_0215_shift1=[select_0215_shift1;select_0215_shift(n,:)]; 
    end
end
figure
    scatter3(select_0215_shift1(:,1), select_0215_shift1(:,2),select_0215_shift1(:,3),1000,select_0215_shift1(:,4),'.');

Color=[1,1,0;
    1,0.666666666666667,0;
    1,0.333333333333333,0;
    1,0,0];
colormap=(sort(Color, 'ascend'));
%caxis([20 80])
h=colorbar;
set(get(h,'label'),'string','Resynchronization time');
xlabel('v_l');ylabel('v_c_2');zlabel('v_c_o_e')
xlim([0 1]);ylim([0.5 2]);zlim([0.5 1]);
set(gca,'fontsize', 18);
%% Plot phase delay index
load select_0215_combineddd.txt
select=[];
for i=1:length(select_0215_combineddd(:,1))
    if select_0215_combineddd(i,6)>0
        select=[select;select_0215_combineddd(i,:)];
    end
end
figure
scatter3(select(:,1), select(:,2),select(:,3),2000,select(:,6),'.')
Color=[1,1,0;
    1,0.666666666666667,0;
    1,0.333333333333333,0;
    1,0,0];
colormap=(sort(Color, 'descend'));
h=colorbar;
%caxis([1 10])
set(get(h,'label'),'string','Phase delay index');
xlabel('v_l');ylabel('v_c_2');zlabel('v_c_o_e')
xlim([0 1]);ylim([0.5 2]);zlim([0.5 1]);
set(gca,'fontsize', 18);
%% Plot alternating shift work
load select_0215_combineddd.txt
select1=[]; select2=[];
for i=1:length(select_0215_combineddd(:,1))
    if select_0215_combineddd(i,4)<=100
    if abs(select_0215_combineddd(i,5))<=0.1
        select1=[select1;select_0215_combineddd(i,:)];
    else
        select2=[select2;select_0215_combineddd(i,:)];
    end
    end
end

figure
scatter3(select1(:,1), select1(:,2),select1(:,3),100,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
hold on
scatter3(select2(:,1), select2(:,2),select2(:,3),100,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
xlabel('v_l');ylabel('v_c_2');zlabel('v_c_o_e')
xlim([0 1]);ylim([0.5 2]);zlim([0.5 1]);
set(gca,'fontsize', 18);
%% Plot transient shift work
load select_0215_combineddd.txt

figure
select=[];
for i=1:length(select_0215_combineddd(:,1))
    if select_0215_combineddd(i,7)<30 & select_0215_combineddd(i,7)>1
        select=[select;select_0215_combineddd(i,:)];
    end
end
scatter3(select(:,1), select(:,2),select(:,3),1000,select(:,7),'.')
Color=[1,1,0;
    1,0.666666666666667,0;
    1,0.333333333333333,0;
    1,0,0];
colormap=(sort(Color, 'descend'));
h=colorbar;
caxis([0 20])
set(get(h,'label'),'string','Maximum phase difference');
xlim([0 1]);ylim([0.5 2]);zlim([0.5 1]);
xlabel('v_l');ylabel('v_c_2');zlabel('v_c_o_e')
set(gca,'fontsize', 18);
%% Correlation between jet lag/shift work properties
load select_0215_combineddd.txt

figure
ttl={'Main view','Left view','Above view','Three-dimensional view'};
select1=[]; select2=[];
for i=1:length(select_0215_combineddd(:,1))
    if select_0215_combineddd(i,4)>1 & select_0215_combineddd(i,4)<100
        if select_0215_combineddd(i,7)>1 
        if abs(select_0215_combineddd(i,5))<0.05 
        select1=[select1;select_0215_combineddd(i,:)];
     
        end
        end
    end
end
for i=1:length(select_0215_combineddd(:,1))
    if select_0215_combineddd(i,4)>1 & select_0215_combineddd(i,4)<100
        if select_0215_combineddd(i,7)>1 
        if abs(select_0215_combineddd(i,5))>0.05 
        select2=[select2;select_0215_combineddd(i,:)];
     
        end
        end
    end
end
plot3(select1(:,4),select1(:,7), select1(:,6),'ro','linewidth',2)
hold on
scatter3(select2(:,4),select2(:,7), select2(:,6),'ko','linewidth',2)
lw=4;
scatter3(select1(:,4),select1(:,7), zeros(length(select1(:,6)),1),'r.','linewidth',lw,'MarkerEdgeAlpha',0.6)
scatter3(select2(:,4),select2(:,7), zeros(length(select2(:,6)),1),'k.','linewidth',lw,'MarkerEdgeAlpha',0.6)
grid on
%xlabel('Resynchronization time (days)')
%ylabel('Maximum phase difference (h)')
%zlabel('Phase delay index (h)')
set(gca,'fontsize', 16);

figure
angle={[0,0],[-90,0],[0 90],[-37.5,30]};
for i=1:4
subplot(2,2,i);

scatter3(select1(:,4), select1(:,7),select1(:,6),200,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
hold on
scatter3(select2(:,4), select2(:,7),select2(:,6),200,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
xlabel('Resynchronization Time (days)')
ylabel('Maximum phase difference (h)')
zlabel('Phase delay index (h)')

xx=0:1:80;yy=0:1:25;zz=0:1:11;
hold on
p=polyfit(select_0215_combineddd(:,4),select_0215_combineddd(:,7),1);
yy1=p(1)*xx+p(2);
plot3(xx,yy1,zeros(1,length(xx)),'b--','linewidth',4)

p=polyfit(select_0215_combineddd(:,7),select_0215_combineddd(:,6),1);
zz1=p(1)*yy+p(2);
plot3(zeros(1,length(yy)),yy,zz1,'b--','linewidth',4)

hold on
p=polyfit(select_0215_combineddd(:,4),select_0215_combineddd(:,6),1);
zz1=p(1)*xx+p(2);
plot3(xx,zeros(1,length(xx)),zz1,'b--','linewidth',4)
x=select_0215_combineddd(:,4);
y=select_0215_combineddd(:,7);
z=select_0215_combineddd(:,6);
lineData=[x,y,z]; 

xyz0(1)=mean(x);
xyz0(2)=mean(y);
xyz0(3)=mean(z); %Fit point coordinates
% Singular value decomposition calculation direction vector
% Singular transformation of covariance matrix
% The direction of the resulting straight line is actually the same as the singular vector corresponding to the largest singular value
 centeredLine=bsxfun(@minus,lineData,xyz0);
 [U,S,V]=svd(centeredLine);
 direction=V(:,1);%Direction vector


t=-100:0.1:100;
xx=xyz0(1)+direction(1)*t;
yy=xyz0(2)+direction(2)*t;
zz=xyz0(3)+direction(3)*t;
xlim([0 80]);ylim([0 25]);zlim([0 10]);
set(gca,'fontsize', 16);
view(angle{i});title(ttl{i});
end
