xcoord = 2.077;

A = [-xcoord .5];
B = [xcoord 1.3];

C = [0.077 -2];
theta = -atan((A(1,2)-B(1,2))/(A(1,1)-B(1,1)));
Rot = [cos(theta) -sin(theta);sin(theta) cos(theta) ];

flat =1;
imp = .5;
if flat ==1
imp = .5;
else
imp = sqrt(2)/2*imp;
end
ycor = -100;

if flat==1
SquareA = [A(1,1)-imp A(1,2)+imp; A(1,1)+imp A(1,2)+imp; A(1,1)+imp A(1,2)-imp; A(1,1)-imp A(1,2)-imp];
SquareB = [B(1,1)-imp B(1,2)+imp; B(1,1)+imp B(1,2)+imp; B(1,1)+imp B(1,2)-imp; B(1,1)-imp B(1,2)-imp];
SquareC = [C(1,1)-imp C(1,2)+imp; C(1,1)+imp C(1,2)+imp; C(1,1)+imp C(1,2)-imp; C(1,1)-imp C(1,2)-imp];
SquareAt = (Rot*SquareA')';
SquareBt = (Rot*SquareB')';
SquareCt = (Rot*SquareC')';
limit = min(SquareAt(:,2));
else 
SquareA = [A(1,1)-imp A(1,2); A(1,1) A(1,2)+imp; A(1,1)+imp A(1,2); A(1,1) A(1,2)-imp];
SquareB = [B(1,1)-imp B(1,2); B(1,1) B(1,2)+imp; B(1,1)+imp B(1,2); B(1,1) B(1,2)-imp];
SquareC = [C(1,1)-imp C(1,2); C(1,1) C(1,2)+imp; C(1,1)+imp C(1,2); C(1,1) C(1,2)-imp];
x = [SquareA(3,1) SquareC(3,1)];
y = [SquareA(3,2) SquareC(3,2)];
mAC = polyfit(x,y,1);
m1 = mAC(1,1);
b1 = mAC(1,2);
x=[SquareB(1,1) SquareC(1,1)];
y = [SquareB(1,2) SquareC(1,2)];
mBC = polyfit(x,y,1);
m2 = mBC(1,1);
b2 = mBC(1,2)
end
while ycor <limit
%while ycor <0
C(1,2) = C(1,2)+.005;
SquareC = [C(1,1)-imp C(1,2)+imp; C(1,1)+imp C(1,2)+imp; C(1,1)+imp C(1,2)-imp; C(1,1)-imp C(1,2)-imp];
SquareCt = (Rot*SquareC')';
x = [SquareAt(2,1) SquareCt(2,1)];
y = [SquareAt(2,2) SquareCt(2,2)];
mAC = polyfit(x,y,1);
m1 = mAC(1,1);
b1 = mAC(1,2);
x = [SquareBt(1,1) SquareCt(1,1)];
y = [SquareBt(1,2) SquareCt(1,2)];
mBC = polyfit(x,y,1);
m2 = mBC(1,1);
b2 = mBC(1,2);
x = [A(1,1),B(1,1)];
y = [A(1,2),B(1,2)];
mAB = polyfit(x,y,1);
m3 = mAB(1,1);
b3 = mAB(1,2);
height = abs(C(1,1)*m3-C(1,2)+b3)/sqrt(m3^2+1);

M = [-m1 1; -m2 1];
D = [b1;b2];
p = (M\D)';
ycor = p(1,2);
%if C(1,2)+.55 >= 0
 %   break;
%end
end
polycoord= [A;B;C];
triangle = polyshape(polycoord);
Area = area(triangle);
%Area = (B(1,1)-A(1,1)*abs(C(1,2)))/2;
distAB = sqrt((A-B)*(A-B)');
distAC = sqrt((A-C)*(A-C)');
distBC= sqrt((C-B)*(C-B)')
avgdistance = (sqrt((A-B)*(A-B)')+sqrt((A-C)*(A-C)')+sqrt((B-C)*(B-C)'))/3;
plot(SquareAt(:,1),SquareAt(:,2)); hold on; plot(SquareBt(:,1),SquareBt(:,2)); hold on; plot(SquareCt(:,1),SquareCt(:,2));
hold on;
fill(SquareA(:,1),SquareA(:,2),'r'); hold on; fill(SquareB(:,1),SquareB(:,2),'b'); hold on; fill(SquareC(:,1),SquareC(:,2),'g');
x = [SquareA(3,1) SquareB(3,1)];
y = [SquareA(3,2) SquareB(3,2)];
plot(x,y);hold on;
x = [SquareA(2,1) SquareC(2,1)];
y = [SquareA(2,2) SquareC(2,2)];
plot(x,y); hold on;
x = [SquareB(1,1) SquareC(1,1)];
y = [SquareB(1,2) SquareC(1,2)];
plot(x,y);
SquareAt = [SquareAt(1,:)-A;SquareAt(2,:)-A;SquareAt(3,:)-A;SquareAt(4,:)-A];
SquareBt = [SquareBt(1,:)-A(1,:); SquareBt(2,:)-A(1,:); SquareBt(3,:)-A;SquareBt(4,:)-A]
%SquareCt = [Ct(1,:)-A(1,:);Ct(2,:)-A(1,:);Ct(3,:)-A;Ct(4,:)-A]
x1 = SquareBt(1,1);
x2 = SquareAt(2,1);
y1 = SquareBt(1,2);
y2 = min(SquareAt(:,2))
b = x2-x1;
c = ((x1-x2)*y1+(y2-y1)*x1)/b;
x1 = [SquareB(1,1) SquareC(1,1)] 
x2 = [SquareA(2,1) SquareC(2,1)];
y1 = [SquareB(1,2) SquareC(1,2)] 
y2 = [SquareA(2,2) SquareC(2,2)];
[xi, yi] = polyxpoly(x1,y1,x2,y2);
scatter(xi,yi);
v_1 = [x1(1,1),x1(1,2),0] - [xi,yi,0];
v_2 = [SquareB(3,1),SquareB(3,2),0] - [xi,yi,0];
ThetaB = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2));
v_1 = [x2(1,1),x2(1,2),0] - [xi,yi,0];
v_2 = [SquareA(3,1),SquareA(3,2),0] - [xi,yi,0];
ThetaA = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2));