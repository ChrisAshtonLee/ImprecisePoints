

%p = [a b c]'


%eq1=p'*inv(A)*p
%eq2 = p'*inv(B)*p
%deltaA = 1;
deltaAp = 1.8;%1.2105;
deltaBp =4;
deltaCp = 1.5;

Ap = [-2 2];
Bp = [6 2];
Cp = [0 -3.5]

altCp = getD(Cp,Ap,Bp);
altBp = getD(Bp,Ap,Cp);
altAp = getD(Ap,Bp,Cp);

C = Cp+[0 0];A = Ap+[0 0]; B = Bp;
altC = getD(C,A,B);
altB = getD(B,A,C);
altA = getD(A,B,C);
deltaA = altA/altAp*deltaAp;
deltaB = altB/altBp*deltaBp;
deltaC = altC/altCp*deltaCp;
theta = 0:pi/16:2*pi;
xa = deltaA*cos(theta)+repmat(A(1),size(theta));
ya = deltaA*sin(theta)+repmat(A(2),size(theta));
xb = deltaB*cos(theta)+repmat(B(1),size(theta));
yb = deltaB*sin(theta)+repmat(B(2),size(theta));
xc = deltaC*cos(theta)+repmat(C(1),size(theta));
yc = deltaC*sin(theta)+repmat(C(2),size(theta));
ratSumOrig =deltaCp/altCp+deltaBp/altBp+deltaAp/altAp;
ratSum = deltaC/altC+deltaB/altB+deltaA/altA;
linesAB = tangentpoints(A,B,deltaA,deltaB);
linesBC = tangentpoints(B,C,deltaB,deltaC);
linesAC = tangentpoints(A,C,deltaA,deltaC);
ABm = linesAB(2,1);  ABb = linesAB(2,2);
ACm = linesAC(1,1); ACb = linesAC(1,2);
BCm = linesBC(2,1);BCb = linesBC(2,2);
ABperp = [1 -1/ABm]/norm([1 -1/ABm]);
ACperp = [1 -1/ACm]/norm([1 -1/ACm]);
BCperp = [1 -1/BCm]/norm([1 -1/BCm]);


ymeet = (ACb-ACm*ABb/ABm)/(1-ACm/ABm);
xmeet = (ymeet-ACb)/ACm;
a = [xmeet ymeet];
ymeet = (ABb-ABm*BCb/BCm)/(1-ABm/BCm);
xmeet= (ymeet-BCb)/BCm;
b = [xmeet ymeet];
ymeet = (ACb-ACm*BCb/BCm)/(1-ACm/BCm);
xmeet = (ymeet-BCb)/BCm;
c = [xmeet ymeet];


vec = [-1 -(ymeet-A(2))/(xmeet-A(1))];
newcenter = A+ 2*vec;
d = abs(newcenter(2)-ACm*newcenter(1)-ACb)/sqrt(1+ACm^2);
xnew = d*cos(theta)+repmat(newcenter(1),size(theta));
ynew = d*sin(theta)+repmat(newcenter(2),size(theta));
altAmod = getD(newcenter,B,C);
newRatio = d/altAmod+deltaB/getD(B,newcenter,C)+deltaC/getD(C,newcenter,B);
figure;
fill(xa,ya,'b','FaceAlpha',0.2)
hold on;
plot([newcenter(1) A(1)],[newcenter(2) A(2)]);
fill(xnew,ynew,'m','FaceAlpha',0.2);
fill(xb,yb,'g','FaceAlpha',0.2)
fill(xc,yc,'r','FaceAlpha',0.2)
scatter([A(1) B(1) C(1)], [A(2) B(2) C(2)],'filled')
ABpm = (A(2)-B(2))/(A(1)-B(1));
ABpb = (A(2)- A(1)*ABpm);
ACpm = (A(2)-C(2))/(A(1)-C(1));
ACpb = A(2)-ACpm*A(1); 
BCpm = (B(2)-C(2))/(B(1)-C(1));
BCpb = B(2) - BCpm*B(1); 

xlim([-6 12]);
ylim([-6 12]);
xplot = -8:.5:12;
yplot = linesAB(2,1)*xplot+linesAB(2,2);
plot(xplot,yplot);
yplot = ABpm*xplot+ABpb;
plot(xplot,yplot);
perpplot = [xplot' yplot']+repmat(-ABperp*deltaC,size(xplot,2),1);
plot(perpplot(:,1),perpplot(:,2),'--');
yplot = ACpm*xplot+ACpb;
plot(xplot,yplot);
perpplot = [xplot' yplot']+repmat(ACperp*deltaB,size(xplot,2),1);
plot(perpplot(:,1),perpplot(:,2),'--');
yplot = linesBC(2,1)*xplot+linesBC(2,2);
plot(xplot,yplot);
yplot = BCpm*xplot+BCpb;
plot(xplot,yplot)
perpplot = [xplot' yplot']+repmat(-BCperp*deltaA,size(xplot,2),1);
plot(perpplot(:,1),perpplot(:,2),'--');

yplot = linesAC(1,1)*xplot+linesAC(1,2);
plot(xplot,yplot);
scatter([a(1) b(1) c(1)],[a(2) b(2) c(2)],'x');
AAArea = .5*abs(a(1)*(b(2)-c(2))+b(1)*(c(2)-a(2))+ c(1)*(a(2)-b(2)));
ocenter = orthocenter([A 0],[B 0],[C 0]);
ocenter = [ocenter(1) ocenter(2)];
odistance_AB= getD(ocenter,A,B);
odistance_AC = getD(ocenter,A,C);
odistance_BC = getD(ocenter,B,C);
scatter(ocenter(1),ocenter(2),'black','filled');
centroid = (a+b+c)/3;
cdistance_AB = getD(centroid,A,B);
cdistance_AC = getD(centroid,A,C);
cdistance_BC = getD(centroid,B,C);

function d = getD(v,P1,P2)
m = (P2(2)-P1(2))/(P2(1)-P1(1));
c = m*P2(1)-P2(2);
a = 1; b = -m; 
d = abs(a*v(2)+b*v(1)+c)/sqrt(a^2+b^2);
end
function lines = tangentpoints(A,B,r1,r2)

a1 = A(1); b1 = A(2); a2 = B(1); b2 = B(2);
  a21 = a2-a1; b21 = b2-b1;
   d2 = a21^2+b21^2;
   r21 = (r2-r1)/d2;
   s21 = sqrt(d2-(r2-r1)^2)/d2; % <-- If d2<(r2-r1)^2, no solution is possible
   u1 = [-a21*r21-b21*s21,-b21*r21+a21*s21]; % Left unit vector
   u2 = [-a21*r21+b21*s21,-b21*r21-a21*s21]; % Right unit vector
   L1 = [a1,b1]+r1*u1; L2 = [a2,b2]+r2*u1; % Left line tangency points
   R1 = [a1,b1]+r1*u2; R2 = [a2,b2]+r2*u2; % Right line tangency points
   points = [L1 L2; R1 R2];
   Lm = (L2(2)-L1(2))/(L2(1)-L1(1));
   Lb = L2(2)-L2(1)*Lm;
   Rm = (R2(2)-R1(2))/(R2(1)-R1(1));
   Rb = R2(2)-R2(1)*Rm;
   lines = [Lm Lb; Rm Rb];
end
function pos = orthocenter(Pa,Pb,Pc)
Pa = Pa(:); Pb = Pb(:); Pc = Pc(:); % Converting to column vectors (if needed)
AB = Pb - Pa; AC = Pc - Pa; BC = Pc - Pb; % Side vectors
N = cross(AC,AB);
L1 = cross(N,BC); L2 = cross(AC,N); % directors
P21 = Pb - Pa;      
P1 = Pa;    
ML = [L1 -L2]; % Coefficient Matrix
lambda = ML\P21;  % Solving the linear system
pos = P1 + lambda(1)*L1; % Line Equation evaluated at lambda(1)
end