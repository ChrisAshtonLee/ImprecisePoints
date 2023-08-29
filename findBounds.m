
function [p1,p2,p3] =  findBounds(hull1,hull2,hull3)
%find bound for hull1 and hull2
constraint12 = zeros(1,3); constraint13 = zeros(1,3); constraint23 = zeros(1,3);
x1 = mean(hull1(:,1)); x2 = mean(hull2(:,1));
y1 = mean(hull1(:,2)); y2 = mean(hull2(:,2));
x3 = mean(hull3(:,1)); y3 = mean(hull3(:,2));
if x1 == x2
    rel12 =[1 2];
end
if y1 == y2 
    rel12 = [1 4];
end
if (x2>x1 && y2>y1) || (x2<x1 && y2<y1)
    rel12 = [1 3];
end
if (x2>x1 && y2<y1) || (x2<x1 && y2>y1)
    rel12 = [2 4];
end
if x2 == x3
    rel23 =[1 2];
end
if y2 == y3 
    rel23 = [1 4];
end
if (x3>x2 && y3>y2) || (x3<x2 && y3<y2)
    rel23 = [1 3];
end
if (x3>x2 && y3<y2) || (x3<x2 && y3>y2)
    rel23 = [2 4];
end
if x1 == x3
    rel13 =[1 2];
end
if y1 == y3 
    rel13 = [1 4];
end
if (x3>x1 && y3>y1) || (x3<x1 && y3<y1)
    rel13 = [1 3];
end
if (x3>x1 && y3<y1) || (x3<x1 && y3>y1)
    rel13 = [2 4];
end
for i =1:4
    if any(rel12(1,:) == i)
    [constraint12(1,1),constraint12(1,2),constraint12(1,3)]= testcorners(hull1(i,:),hull2(i,:),hull3,i);
    end
    if ~isequal(constraint12,[0 0 0])
        break;
    end
end
for i =1:4
    if any(rel13(1,:) == i)
    [constraint13(1,1),constraint13(1,2),constraint13(1,3)]= testcorners(hull1(i,:),hull3(i,:),hull2,i);
    end
    if ~isequal(constraint13,[0 0 0])
        break;
    end
end
for i =1:4
    if  any(rel23(1,:) == i)
    [constraint23(1,1),constraint23(1,2),constraint23(1,3)]= testcorners(hull2(i,:),hull3(i,:),hull1,i);
    end
    if ~isequal(constraint23,[0 0 0])
        break;
    end
end
%check if feasible
if isequal(constraint12,[0 0 0]) || isequal(constraint13,[0 0 0])|| isequal(constraint23,[0 0 0])
    p1 = 0; p2 = 0; p3 = 0;
    
else
    A = [-constraint12(1,1:2);-constraint13(1,1:2);-constraint23(1,1:2)];
    B = [-constraint12(1,3);-constraint13(1,3);-constraint23(1,3)];
    f = ones(1,2);
    %if any(logical(linprog(f,A,B)))
        A1213= [constraint12(1,1:2); constraint13(1,1:2)];
        B1213= [constraint12(1,3);constraint13(1,3)];
        p1 = fliplr((inv(A1213)*B1213)');
        %p1 = fliplr(linsolve(A1213,B1213)');
        A1323= [constraint13(1,1:2); constraint23(1,1:2)];
        B1323= [constraint13(1,3);constraint23(1,3)];
        p2 = fliplr((inv(A1323)*B1323)');
        %p2 = fliplr(linsolve(A1323,B1323)');
        A1223= [constraint12(1,1:2); constraint23(1,1:2)];
        B1223= [constraint12(1,3);constraint23(1,3)];
        p3 = fliplr((inv(A1223)*B1223)');
        %p3 = fliplr(linsolve(A1223,B1223)');
        testp = [(p1(1,2)+p2(1,2)+p3(1,2))/3; (p1(1,1)+p2(1,1)+p3(1,1))/3];
        if any(A*testp > B)
   
        p1 = 0; p2 = 0; p3 = 0;
        end
end


function [a1,a2,a3] = testcorners(hullA,hullB,hullC,type)
a1 = 0;a2=0;a3 = 0;
switch type
    case 1 
        xdir = -1; ydir = 1; xydir = 1;
    case 2
        xdir = 1; ydir = 1; xydir = 1;
    case 3
        xdir = 1; ydir = -1; xydir = -1;
    case 4
        xdir = -1; ydir = -1; xydir = -1;
end
x = [hullA(1,1) hullB(1,1)];
y = [hullA(1,2) hullB(1,2)];
%check xbound
if hullA(1,1) == hullB(1,1)
  x =hullA(1,1);
  if (xbound(x,hullC,xdir))
      a1 = 0;a2 = xdir; a3 = xdir*x;
  end
  %check ybound
elseif hullA(1,2) == hullB(1,2)
    y = hullA(1,2);
      if (ybound(y,hullC,ydir))
      a1 = ydir; a2 = 0; a3 = ydir*y;
      end
else
p = polyfit(x,y,1);
m = p(1,1);
b = p(1,2);
    if xybound(m,b,hullC,xydir)
        a1 = xydir; a2 = -xydir*m; a3 = xydir*b;
    end
end

end
function innerbound = xbound(x,hull,dir)
innerbound = 0;
switch dir
    case 1
        if ~any(hull(:,1)<x)
            innerbound = 1;
        end
    case -1
        if ~any(hull(:,1)>x)
            innerbound = 1;
        end
end       
end
function innerbound = ybound(y,hull,dir)
innerbound = 0;
switch dir
    case 1
        if ~any(hull(:,2)<y)
            innerbound = 1;
        end
    case -1
        if ~any(hull(:,2)>y)
            innerbound = 1;
        end
end       
end
function innerbound = xybound(m,b,hull,dir)
switch dir
    case 1
        if ~any(hull(:,2)< m*hull(:,1)+b)
            innerbound = 1;
        else
            innerbound = 0;
        end
    case -1
        if ~any(hull(:,2)> m*hull(:,1)+b)
            innerbound = 1;
        else
            innerbound = 0;
        end        
end
end
end
