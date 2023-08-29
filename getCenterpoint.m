function [c] = getCenterpoint(X)
Xf = X;
n = size(Xf,1);
k = ceil(2/3*n)+1;
comb = combnk(1:n,k);
tuple = zeros(k,2,size(comb,1));

for i = 1:size(comb,1)
    for j = 1:k
        tuple(j,:,i) = Xf(comb(i,j),:);
    end
    
end
try
I = intersectionHull('vert',tuple(:,:,1),'vert',tuple(:,:,2),'vert',tuple(:,:,3),'vert',tuple(:,:,4),'vert',tuple(:,:,5),'vert',tuple(:,:,6));
c = I.vert;
catch
    warning('centerpoint failed');
    c = Xf;
end

end