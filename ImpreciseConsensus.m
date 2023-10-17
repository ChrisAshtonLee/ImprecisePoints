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
dt=0.008; % numerical steplength
Tf=5; % final time
network_no=1;

IHalg = 0;
%imp = 0;
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
%attacker = 10;
iterlimit = 2500;

count = zeros(N,iterlimit);
t=0; 
imp = 1;
for imp = 1:1
[X,n,N]= load_network(network_no,Delta);
num_agents = N; 
attacker = N;
attacker = 10;
init_pos = X;
attacktarget = [8 ;8];
Xhist = zeros(n,num_agents,iterlimit);
safepoint= zeros(2,N);
d = 2;
iter = 1;
dmin = inf;
meanMovement = inf;
Movement = zeros(1,num_agents);

while (meanMovement >.001 && iter <= iterlimit)
DX=zeros(n,N); 

for i=1:num_agents
    Est = X;
    for j = 1:N
        if j ~= i
            Est(:,j) = Est(:,j)+[2*imp*(rand()-.5);2*imp*(rand()-.5)];
            if i<3 & attacker == 6
                Est(1,attacker)  = -X(1,attacker);
            elseif i == 3 & iter<500 & attacker ==6
                Est(1,attacker) = -X(1,attacker);
            end
        end
    end
    [safepoint(:,i)] = BadCenterpoint(X,Est,i);
    %safepoint(:,i) = DiscrepancyPoint(X,Est)';
    if IHalg
    cp_n = ceil(d/(d+1)*num_agents)+1;
    cpEst = zeros(2,cp_n);
    cp_indices = combnk(1:num_agents,cp_n);
    num_cp_vertices = 3*nchoosek(cp_n,3);
    %IH_cpEst = zeros( num_cp_vertices,2,size(cp_indices,1));
    IH_cpEst ={}
    testIH = [];
    for cp_i = 1: size(cp_indices,1)
    cpEst = Est(:,cp_indices(cp_i,:));
    [~,tempcpEst] = getInvariant(cpEst,imp,i);
    dim = size(tempcpEst,1);
    IH_cpEst{cp_i} = tempcpEst(1:dim,:);
    end;
    try
    I = intersectionHull('vert',IH_cpEst{1},'vert',IH_cpEst{2},'vert',IH_cpEst{3},'vert',IH_cpEst{4},'vert',IH_cpEst{5},'vert',IH_cpEst{6});
    CP_vertices= I.vert;
    safepoint(:,i) = [mean(CP_vertices(:,1)); mean(CP_vertices(:,2))];
    catch
        safepoint(:,i) = X(:,i);
    end
    end
    %safepoint(:,i) = getInvariant(Est,imp,i);
end

%% Update the states using an Euler approximation
    for i=1:num_agents
     if i == attacker
       X(:,i) =X(:,i)+ dt*5.*(attacktarget-X(:,i)); 
     else
        DX(:,i)=DX(:,i)+(safepoint(:,i)-X(:,i)); 
        X(:,i)=X(:,i)+dt.*DX(:,i);
        Movement(1,i) =  dot([1 1],dt.*DX(:,i));
     end
    end
 MeanMovement = sum(Movement(1,:))/(num_agents-1);

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
   x = reshape(Xhist(1,:,:),num_agents,iterlimit);
   y = reshape(Xhist(2,:,:),num_agents,iterlimit);
   t = 1:iterlimit;
   xupperB = max(init_pos(1,:)).*ones(1,iterlimit);
   xlowerB = min(init_pos(1,:)).*ones(1,iterlimit);
   yupperB = max(init_pos(2,:)).*ones(1,iterlimit);
   ylowerB = min(init_pos(2,:)).*ones(1,iterlimit);
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
save('IMP= '+string(imp));
end
function [target] = BadCenterpoint(X,Est,agent)
Xf = X';
Estf = Est';
Xavg = [mean(Xf(:,1)) mean(Xf(:,2))];
c = getCenterpoint(Estf,agent);
vertices = convhull(Xf);
truehull = Xf(vertices,:);
%target = c(1,:)';
target = [mean(c(:,1)) mean(c(:,2))];

max = 0;

for i = 1:size(c,1)
   % if ~(inpolygon(c(i,1),c(i,2),truehull(:,1),truehull(:,2)))
        dist = sqrt(abs(dot(c(i,:)',Xavg)));
        if  dist  >max
            max = dist;
            target = c(i,:)';
        end
       
    %end
end

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
