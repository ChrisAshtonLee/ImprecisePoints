x = [init_pos(1,:),init_pos(1,1)];
y = [init_pos(2,:),init_pos(1,2)];

plot(x,y,'b');
axis([-5 15 -10 10])
hold on;
for i = 1:6
    x = reshape(Xhist(1,i,:),1,1000);
    y = reshape(Xhist(2,i,:),1,1000);
   plot(x,y);
   hold on; 
end
for i = 1:6
    x = Xhist(1,i,size(Xhist,3));
    y = Xhist(2,i,size(Xhist,3));
    plot(x,y,'x');
    hold on;
end
legend('conv(init pos)','agent 1','agent 2',' agent 3', 'agent 4', 'agent 5', 'attacker');
xlabel('x pos');
ylabel('y pos');
title('5 agents 1 attacker(non-Byzantine)');
