%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

figure
%% Proximity disk
Delta=3; % interaction distance

% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions

mode = 5;
%% Some numerical integration parameters
dt=0.01; % numerical steplength
Tf=5; % final time
network_no=1;
imp = 1;
plothist =1;
if mode == 1
   network_no=1; 
   [X,n,N]= load_network(network_no,Delta);
   num_agents = N; 
   attacker = N;
   init_pos = X;
end
if mode == 4
   network_no=4; 
   [X,n,N]= load_network(network_no,Delta);
   num_agents = N; 
   attacker = N;
   init_pos = X;
   attacktarget = [8 ;8];
end
if mode == 5
   network_no=5; 
   [X,n,N]= load_network(network_no,Delta);
   num_agents = N; 
   attacker = N;
   init_pos = X;
   attacktarget = [8 ;8];
end
Xhist = zeros(n,num_agents,1000);
safepoint= zeros(2,N);
count = zeros(N,1000);
t=0; 
dist = 0;
iter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step = 0.01;
%dist = @(xi, xj) sqrt( (xi - xj)' * (xi - xj));
%weight = @(xi, xj) weightfcn(dist(xi, xj)); 
% H for hysterisys
epsilon = .05;
dmin = inf;
%% A is the adjacency matrix associated with the system
A = zeros(N)
for i =1:(N-1)
    for j = (i+1):N
        A(i,j) = 1;
    end
end

partition = zeros(N,2);
while (iter<=1000);
DX=zeros(n,N); 

for i=1:num_agents
    Est = X;
    for j = 1:N
        if j ~= i
            Est(:,j) = Est(:,j)+[2*imp*(rand()-.5);2*imp*(rand()-.5)];
        end
    end
    [safepoint(:,i), count(i,iter)] = BadCenterpoint(X,Est);
    %safepoint(:,i) = DiscrepancyPoint(X,Est)';
    
    %safepoint(:,i) = getInvariant(Est,imp);
end

%% Update the states using an Euler approximation
    for i=1:num_agents
     if i == attacker
       X(:,i) =X(:,i)+ dt.*(attacktarget-X(:,i)); 
     else
        DX(:,i)=DX(:,i)+(safepoint(:,i)-X(:,i)); 
        X(:,i)=X(:,i)+dt.*DX(:,i);
     end
    end
 

%% Update time
  t=t+dt;
  Xhist(:,:,iter) = X;

%% Plot the solution every 10 iterations
  if (mod(iter,10)==0);
    plotsol(X,N,imp,init_pos,safepoint);
  end;

  iter=iter+1;
end;
if plothist == 1
   x = reshape(Xhist(1,:,:),num_agents,1000);
   y = reshape(Xhist(2,:,:),num_agents,1000);
   t = 1:1000;
   xupperB = max(init_pos(1,:)).*ones(1,1000);
   xlowerB = min(init_pos(1,:)).*ones(1,1000);
   yupperB = max(init_pos(2,:)).*ones(1,1000);
   ylowerB = min(init_pos(2,:)).*ones(1,1000);
   figure;
   plot(t,x);hold on;
   plot(t,xupperB,'--'); plot(t,xlowerB,'--')
   xlabel('Iterations');
   ylabel('Position(X)');
   
   figure;
   plot(t,y);hold on;
   plot(t,yupperB,'--'); plot(t,ylowerB,'--');
   xlabel('Iterations');
   ylabel('Position(Y)');
end
function [target,count] = BadCenterpoint(X,Est)
Xf = X';
Estf = Est';
Xavg = [mean(Xf(:,1)) mean(Xf(:,2))];
c = getCenterpoint(Estf);
vertices = convhull(Xf);
truehull = Xf(vertices,:);
%target = c(1,:)';
target = [mean(c(:,1)) mean(c(:,2))];
count = 0;
max = 0;
%{
for i = 1:size(c,1)
    if ~(inpolygon(c(i,1),c(i,2),truehull(:,1),truehull(:,2)))
        dist = sqrt(c(i,:)'*Xavg);
        if  dist  >max
            max = dist;
            target = c(i,:)';
        end
        count = count+1;
    end
end
%}
end
function target =  DiscrepancyPoint(X,Est)
Xf = X';
Estf = Est';
vertices = convhull(Xf);
truehull = Xf(vertices,:);
vertices = convhull(Estf);
esthull = Estf(vertices,:);
%esthull = getCenterpoint(Estf);
target = [mean(Estf(:,1)) mean(Estf(:,2))];
Xavg = [mean(Xf(:,1)) mean(Xf(:,2))];
max = 0;

found = 0;
try
I = intersectionHull('vert',truehull,'vert',esthull);
found = size(I.vert,1);
catch
    warning('no intersection');
  
end
if found~= 0
for i = 1:size(esthull,1)-1
    
    if ~any(I.vert == esthull(i,:))
        dist = sqrt(esthull(i,:)'*Xavg);
        if dist > max
            max = dist;
            target = esthull(i,:);
        end
    end
end
end
end
