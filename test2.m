Ax = A(1); Ay= A(2);
Bx = B(1); By = B(2);
Cx = C(1); Cy = C(2);
ra = .5;
rb = 1.5;
rc = 1;
raplus = 1;
rbplus = 1;
rcplus = 1.5;

Mab = (By-Ay-(rb-ra))/(Bx-Ax+(rb-ra));Nab = 1/Mab;
Mac = (Cy-Ay+(rc-ra))/(Cx-Ax+(rc-ra)); Nac = 1/Mac;
Mbc = (By-Cy+(rb-rc))/(Bx-Cx-(rb-rc)); Nbc = 1/Mbc;

%Ax = ra*(2-Mac+Mab)/(Mac-Mab); Ay = -ra*(Nac+Nab)/(Nac-Nab);
%Bx = rb*(2+Mbc+Mab)/(Mbc-Mab); By = rb*(2+Nab+Nbc)/(Nab-Nbc);
%Cx = rc*(Mbc+Mac)/(Mbc-Mac); Cy = rc*(2+Nbc-Nac)/(Nac-Nbc);
Axfac = (2-Mac+Mab)/(Mac-Mab); Ayfac = (Nac+Nab)/(Nac-Nab);
Bxfac = (2+Mbc+Mab)/(Mbc-Mab); Byfac = (2+Nab+Nbc)/(Nab-Nbc);
Cxfac = (Mbc+Mac)/(Mbc-Mac); Cyfac = (2+Nbc-Nac)/(Nac-Nbc);
Axplus = (raplus-ra)*Axfac+Ax; Ayplus =  -(raplus-ra)*Ayfac+Ay;
Bxplus = (rbplus-rb)*Bxfac+Bx; Byplus = (rbplus-rb)*Byfac+By;
Cxplus = (rcplus-rc)*Cxfac+Cx; Cyplus = (rcplus-rc)*Cyfac+Cy;
Asquareplus = makeSquare([Axplus Ayplus],raplus);
Bsquareplus = makeSquare([Bxplus Byplus],rbplus);
Csquareplus = makeSquare([Cxplus Cyplus],rcplus);

scatter([Axplus Bxplus Cxplus],[Ayplus Byplus Cyplus]);
fill(Asquareplus(:,1),Asquareplus(:,2), 'r','FaceAlpha',0.3);
fill(Bsquareplus(:,1),Bsquareplus(:,2), 'b','FaceAlpha',0.3);
fill(Csquareplus(:,1),Csquareplus(:,2), 'g','FaceAlpha',0.3);
function [Square] = makeSquare(A,imp)
Square = [A(1,1)-imp A(1,2)+imp; A(1,1)+imp A(1,2)+imp; A(1,1)+imp A(1,2)-imp; A(1,1)-imp A(1,2)-imp];

end