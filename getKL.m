function [KL,w] = getKL(signal,label,from)
w = zeros(1,size(label,2));
for qidx = 1:size(label,2)
    q = qidx+from-1;
    for r=1:q
        clstr = signal(:,label(:,qidx)==r);
        D(r) = 0;
        n(r) = size(clstr,2);
        clstrsq = dot(clstr,clstr,1);
        D(r) = sum(sum(bsxfun(@plus,clstrsq',clstrsq)-2*(clstr'*clstr)));
    end
    w(qidx) = (1./(2*n)) * D';
    M(qidx) = w(qidx).*q^(2/size(signal,1));
end
d = M(1:end-1)-M(2:end);
KL = (d(1:end-1) - d(2:end))./M(1:end-2);
KL(d(1:end-1)<0) = 0;
KL(d(1:end-1)<d(2:end)) = 0;
end