plotImp = 1;
agent6 = 'agent 6';
%agent6 = 'adversary';
f = figure;

f.Position=[10 10 1000 800];
x = [init_pos(1,1:6),init_pos(1,1)];
y = [init_pos(2,1:6),init_pos(1,2)];
plot(x,y,'b--','LineWidth',1.5);
axis([-5 15 -10 10]);
hold on;



for i = 1:6
    x = reshape(Xhist(1,i,1:iterlimit),1,iterlimit);
    y = reshape(Xhist(2,i,1:iterlimit),1,iterlimit);
   plot(x,y,'LineWidth',1.5);
   hold on; 
end
if plotImp
    for i=1:6
    x = [Xhist(1,i,end)-imp Xhist(1,i,end)+imp Xhist(1,i,end)+imp Xhist(1,i,end)-imp Xhist(1,i,end)-imp]; y =[Xhist(2,i,end)+imp Xhist(2,i,end)+imp Xhist(2,i,end)-imp Xhist(2,i,end)-imp Xhist(2,i,end)+imp]; 
    plot(x,y,'--r', 'LineWidth',1.5);
    end
    hold on;
end
for i = 1:6
  plot(Xhist(1,i,iterlimit),Xhist(2,i,iterlimit),'x','LineWidth',1.5,'MarkerSize',10);
end
%{
x = -1.*reshape(Xhist(1,attacker,:),1,iterlimit);
    y = reshape(Xhist(2,attacker,:),1,iterlimit);
   plot(x,y,'LineWidth',1.5);
   hold on; 
%}

if plotImp
    lgd = legend('conv(X_0)','agent 1','agent 2','agent 3', 'agent 4', 'agent 5', agent6,'potential regions');
    lgd.FontSize = 24;
else
    lgd = legend('conv(X_0)','agent 1','agent 2','agent 3', 'agent 4', 'agent 5', agent6);
    lgd.FontSize = 24;
end
xlabel('x pos');
ylabel('y pos');
%title('5 Agents 1 Adversary with Imprecision');
%title('5 Agents 1 Adversary with Imprecision.');
