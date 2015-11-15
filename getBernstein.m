function [B]=getBernstein(p,xi)

numqpt = length(xi);
B  = zeros(p+1,numqpt);
for i=1:p+1
    n = i-1;
    B(i,:) = nchoosek(p,n)*(1-0.5*(xi+1)).^(p-n).*(0.5*(xi+1)).^n;
end

