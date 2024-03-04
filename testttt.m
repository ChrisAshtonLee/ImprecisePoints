xcoord = 4;

A = [-xcoord .5];
B = [xcoord 2.5];
%acn= ( (A+[1 -1])-(B+[2 -2])+( C+[-1 1])-(B+[-2 2]) )/2;
%translate =acn*sqrt(2)/sqrt(8);

C = [1 -4];
theta = -atan((A(1,2)-B(1,2))/(A(1,1)-B(1,1)));
Rot = [cos(theta) -sin(theta);sin(theta) cos(theta) ];

flat =1;
imp = .5;
if flat ==1
imp = 1;
else
imp = sqrt(2)/2*imp;
end
ycor = -100;

if flat==1
impB = 2;
%SquareA = [A(1,1)-imp A(1,2)+imp; A(1,1)+imp A(1,2)+imp; A(1,1)+imp A(1,2)-imp; A(1,1)-imp A(1,2)-imp];
SquareA = makeSquare(A,.5*imp);
SquareB = makeSquare(B,1.5*imp);
%SquareB = [B(1,1)-impB B(1,2)+impB; B(1,1)+impB B(1,2)+impB; B(1,1)+impB B(1,2)-impB; B(1,1)-impB B(1,2)-impB];
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
b2 = mBC(1,2);

end
%limit = -2.995;
while ycor <limit-.8
%while C(1,2) <-2.1399
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
distBC= sqrt((C-B)*(C-B)');

avgdistance = (sqrt((A-B)*(A-B)')+sqrt((A-C)*(A-C)')+sqrt((B-C)*(B-C)'))/3;
v1 = [A(1,1) A(1,2) 0]-[C(1,1) C(1,2) 0];
v2 = [B(1,1) B(1,2) 0]-[C(1,1) C(1,2) 0];
ThetaACB= atan2d(norm(cross(v1,v2)),dot(v1,v2));
v1 = [C(1,1) C(1,2) 0]-[A(1,1) A(1,2) 0];
v2 = [B(1,1) B(1,2) 0]-[A(1,1) A(1,2) 0];
ThetaBAC= atan2d(norm(cross(v1,v2)),dot(v1,v2));
v1 = [C(1,1) C(1,2) 0]-[B(1,1) B(1,2) 0];
v2 = [A(1,1) A(1,2) 0]-[B(1,1) B(1,2) 0];
ThetaABC= atan2d(norm(cross(v1,v2)),dot(v1,v2));
%plot(SquareAt(:,1),SquareAt(:,2)); hold on; plot(SquareBt(:,1),SquareBt(:,2)); hold on; plot(SquareCt(:,1),SquareCt(:,2));
%hold on;
AlphaA = -atan((C(1,2)-A(1,2))/(C(1,1)-A(1,1)));
BetaA = ThetaBAC/360*2*pi;
CharlieA = pi-AlphaA-BetaA;
ACnew= sqrt(2)*sin(CharlieA)/sin(BetaA);
point = A+(B-A)/norm(B-A)*ACnew;

fill(SquareA(:,1),SquareA(:,2),'r','FaceAlpha',0.6); hold on; fill(SquareB(:,1),SquareB(:,2),'b','FaceAlpha',0.6); hold on; fill(SquareC(:,1),SquareC(:,2),'g','FaceAlpha',0.6);
x = [SquareA(3,1) SquareB(3,1)];
y = [SquareA(3,2) SquareB(3,2)];
plot(x,y);hold on;
x = [SquareA(2,1) SquareC(2,1)];
y = [SquareA(2,2) SquareC(2,2)];
plot(x,y); hold on;
x = [SquareB(1,1) SquareC(1,1)];
y = [SquareB(1,2) SquareC(1,2)];
plot(x,y);
x = [A(1,1) B(1,1) C(1,1) A(1,1)];
y = [A(1,2) B(1,2) C(1,2) A(1,2)];
plot(x,y);

b1= 1;
b2 =1;
BCf = (-(B+[-1.5 1.5]) +(C+[-1 1]));
a1 = -BCf(2)/BCf(1);
n1 = -1/a1;
c1 =  -(B(1,2)+1.5+a1*(B(1,1)-1.5));
BAf = (-(B+[1.5 -1.5]) + (A+[.5 -.5]));
a2 = -BAf(2)/BAf(1);
n2 = -1/a2;
c2 = -(B(1,2)-1.5+a2*(B(1,1)+1.5));
point = [(b1*c2-b2*c1)/(a1*b2-a2*b1) (a2*c1-a1*c2)/(a1*b2-a2*b1)];
vector = point-B;
Bmod = B+2/3*vector;
a1 = -a1;
a2 = -a2;
x1 = (3+a2*3)/(a1-a2);
x0 = x1/2;
y_1 = (3+3*n2)/(n2-n1);
y_0 = (2+2*n2)/(n2-n1);
y_1 = (3/(n2-n2^2+n1*n2-n1))

x1 = (3+1.5*a1+1.5*a2)/(a1-a2);
x0 = x1*2/3;
y_1 = (3+1.5*n2+1.5*n1)/(n2-n1);
y_0 = y_1*2/3;
answer  = x1-x0;
answer2 = y_1-y_0;
BCf = (-(B+[-1.5 1.5]) +(C+[-1 1]));
Bmod = B-[answer answer2];
ACf = A+[.5 .5]-C-[1 1];
a1 = ACf(2)/ACf(1);
n1 = 1/a1;
x1 = (1-.5*a1+.5*a2)/(a1-a2);
x0 = x1*2;
y_1 = (.5*n2+.5*n1)/(n2-n1);
y_0 = y_1*2;
ccircrad= distAB*distBC*distAC/sqrt((distAB+distAC+distBC)*(-distAB+distAC+distBC)*(distAB-distAC+distBC)*(distAB+distAC-distBC));
answer  = x1-x0;
answer2 = y_1-y_0;
Amod = A-[answer answer2];
Amodsquare = makeSquare(Amod,1);
fill(Amodsquare(:,1),Amodsquare(:,2),'r','FaceAlpha',0.3);
c1 =  -(B(1,2)+2+a1*(B(1,1)-2));
BAf = (-(B+[1.5 -1.5]) + (A+[.5 -.5]));
Bmodsquare = makeSquare(Bmod,1);
fill(Bmodsquare(:,1),Bmodsquare(:,2), 'b','FaceAlpha',0.3)
%{
SquareA4 = makeSquare(SquareA(4,:),imp);
SquareB4 = makeSquare(SquareB(4,:),imp);
SquareC4 = makeSquare(SquareC(4,:),imp);
fill(SquareA4(:,1),SquareA4(:,2),'r','FaceAlpha',0.3); hold on; fill(SquareB4(:,1),SquareB4(:,2),'b','FaceAlpha',0.3); hold on; fill(SquareC4(:,1),SquareC4(:,2),'g','FaceAlpha',0.3);
[x4,y4] = getIH(SquareA4,SquareB4,SquareC4);
scatter(x4,y4);hold on;
SquareA1 = makeSquare(SquareA(1,:),imp);
SquareB1 = makeSquare(SquareB(1,:),imp);
SquareC1 = makeSquare(SquareC(1,:),imp);
fill(SquareA1(:,1),SquareA1(:,2),'r','FaceAlpha',0.3); hold on; fill(SquareB1(:,1),SquareB1(:,2),'b','FaceAlpha',0.3); hold on; fill(SquareC1(:,1),SquareC1(:,2),'g','FaceAlpha',0.3);
[x1,y1] = getIH(SquareA1,SquareB1,SquareC1);
scatter(x1,y1);hold on;
SquareA2 = makeSquare(SquareA(2,:),imp);
SquareB2 = makeSquare(SquareB(2,:),imp);
SquareC2 = makeSquare(SquareC(2,:),imp);
fill(SquareA2(:,1),SquareA2(:,2),'r','FaceAlpha',0.3); hold on; fill(SquareB2(:,1),SquareB2(:,2),'b','FaceAlpha',0.3); hold on; fill(SquareC2(:,1),SquareC2(:,2),'g','FaceAlpha',0.3);
[x2,y2] = getIH(SquareA2,SquareB2,SquareC2);
scatter(x2,y2);hold on;
%}
%{
SquareAf = makeSquare(SquareA(1,:),imp);
SquareBf = makeSquare(SquareB(2,:),imp);
for theta = pi:pi/8:2*pi
SquareCf = makeSquare(C+1*[cos(theta) sin(theta)],imp);
%fill(SquareAf(:,1),SquareAf(:,2),'r','FaceAlpha',0.3); hold on; fill(SquareBf(:,1),SquareBf(:,2),'b','FaceAlpha',0.3); hold on; fill(SquareCf(:,1),SquareCf(:,2),'g','FaceAlpha',0.3);
fill(SquareCf(:,1),SquareCf(:,2),'g', 'FaceAlpha',0.3);
[xf,yf] = getIH(SquareA,SquareB,SquareCf);
scatter(xf,yf);hold on;
end
%}
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
x1 = [SquareA(2,1) SquareC(2,1)] 
x2 = [A(1,1) B(1,1)];
y1 = [SquareA(2,2) SquareC(2,2)] 
y2 = [A(1,2) B(1,2)];
[xac, yac] = polyxpoly(x1,y1,x2,y2);
x1 = [SquareB(1,1) SquareC(1,1)] 
x2 = [A(1,1) B(1,1)];
y1 = [SquareB(1,2) SquareC(1,2)] 
y2 = [A(1,2) B(1,2)];
[xbc, ybc] = polyxpoly(x1,y1,x2,y2);
sqrt((A-B)*(A-B)');
dsim = sqrt(([xac yac]-[xbc ybc])*([xac yac]-[xbc ybc])');
leftover = distAB-2*dsim;
scatter([xi xac xbc],[yi yac ybc]);
polycoord= [[xi yi];[xac yac];[xbc ybc]];
triangle = polyshape(polycoord);
Area2 = area(triangle);
AreaRatio = Area2/Area;
AB_modlength = sqrt((xac-xbc)^2+(yac-ybc)^2);
ABbyDAB = distAB/AB_modlength;
v_1 = [x1(1,1),x1(1,2),0] - [xi,yi,0];
v_2 = [SquareB(3,1),SquareB(3,2),0] - [xi,yi,0];
ThetaB = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2));
v_1 = [x2(1,1),x2(1,2),0] - [xi,yi,0];
v_2 = [SquareA(3,1),SquareA(3,2),0] - [xi,yi,0];
ThetaA = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2));
function [xi,yi]= getIH(SquareA,SquareB,SquareC)
x1 = [SquareB(1,1) SquareC(1,1)] 
x2 = [SquareA(2,1) SquareC(2,1)];
y1 = [SquareB(1,2) SquareC(1,2)] 
y2 = [SquareA(2,2) SquareC(2,2)];
[xi, yi] = polyxpoly(x1,y1,x2,y2);
end
function [Square] = makeSquare(A,imp)
Square = [A(1,1)-imp A(1,2)+imp; A(1,1)+imp A(1,2)+imp; A(1,1)+imp A(1,2)-imp; A(1,1)-imp A(1,2)-imp];

end