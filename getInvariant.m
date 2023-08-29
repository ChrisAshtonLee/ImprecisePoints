function [IH] = getInvariant(X,imp)
Xf = X';
n = size(Xf,1);
Xfhull = zeros(4,2,n);
for i = 1:n
 Xfhull(1,:,i) =[Xf(i,1)-imp Xf(i,2)+imp];
 Xfhull(2,:,i) = [Xf(i,1)+imp Xf(i,2)+imp];
 Xfhull(3,:,i) = [Xf(i,1)+imp Xf(i,2)-imp];
 Xfhull(4,:,i) = [Xf(i,1)-imp Xf(i,2)-imp];
end
triples = combnk(1:n,3);
num_triples = size(triples,1);
IHtriples = zeros(3*num_triples,2);
for i= 1:num_triples
   hull1 = Xfhull(:,:,triples(i,1));
   hull2 = Xfhull(:,:,triples(i,2));
   hull3 = Xfhull(:,:,triples(i,3));
   [IHtriples(3*i-2,:), IHtriples(3*i-1,:),IHtriples(3*i,:)] = findBounds(hull1,hull2,hull3);
end
try
IHindices = convhull(IHtriples);
IH = IHtriples(IHindices,:);
IH = [mean(IH(:,1)) mean(IH(:,2))];
catch
    warning('IH not found');
    IH = [mean(Xf(:,1)) mean(Xf(:,2))];
end
end
