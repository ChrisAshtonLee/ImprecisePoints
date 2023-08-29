function[]=show_graph(X,N,imp,init_pos,safepoint);

%figure(1);

Delta = 3;
hold off;
for i =1:N;
plot(init_pos(1,i),init_pos(2,i),'diamond','HandleVisibility','off');
hold on;
end
x = [init_pos(1,:),init_pos(1,1)];
y = [init_pos(2,:),init_pos(1,2)];
plot(x,y,'b');
hold on;
for i=1:N;
  x = [X(1,i)-imp X(1,i)+imp X(1,i)+imp X(1,i)-imp X(1,i)-imp]; y =[X(2,i)+imp X(2,i)+imp X(2,i)-imp X(2,i)-imp X(2,i)+imp]; 
plot(x,y,'r', 'LineWidth',1);
hold on;
  plot(X(1,i),X(2,i),'o');
  hold on;
  plot(safepoint(1,i),safepoint(2,i),'x');
  hold on;
  plot([X(1,i),safepoint(1,i)],[X(2,i),safepoint(2,i)],'HandleVisibility','off');
end; 


axis([-4*Delta,4*Delta,-4*Delta,4*Delta]);
drawnow;

