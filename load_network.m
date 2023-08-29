function[X,n,N]=load_network(network_no,Delta);
n=2; X=[];
if (network_no==1);
  N=4;  
  Delta=3;
  X(:,1)=[-Delta;-Delta]; X(:,2)=[-Delta;Delta];
  X(:,3)=[Delta;Delta]; X(:,4)=[Delta,-Delta];

elseif (network_no==2);
  N=8;
  Delta = 3;
  for i = 1:N-2
     X(:,i) = Delta*cos(pi([-10 10],2,1));
  end
  X(:,N-1) = [-4;6];
  X(:,N) = [4;6];
elseif (network_no == 4);
  N = 7;
  X(:,1) = [0;0];
  X(:,2) = [0;1];
  X(:,3) = [1;0];
  X(:,4) = [5;2];
  X(:,5) = [-2;1];
  X(:,6) = [2;-2];
  X(:,7) = [-2;6];
elseif(network_no == 5);
    N = 6;
    X(:,1) = [-4;-4];
    X(:,2) = [-4;0];
    X(:,3) = [-3;3];
    X(:,4) = [1;4];
    X(:,5) = [4;0];
    X(:,6) = [4;-3];
else
  N=7;
  Delta=1*Delta;
  X(:,1)=[-8*Delta/5;0]; X(:,2)=[-4*Delta/5;0];
  X(:,3)=[0;0]; X(:,4)=[0;4*Delta/5]; X(:,5)=[0;8*Delta/5];
  X(:,6)=[4*Delta/5;0]; X(:,7)=[8*Delta/5;0];
end;
