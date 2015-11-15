function [dB]=getdBernstein(p,xi)

numqpt = length(xi);
dB  = zeros(p+1,numqpt);
B = getBernstein(p-1,xi);
for i=1:p+1
    n = i-1;
    if(n == 0)
        dB(i,:) = 0.5*p*(-B(i,:));
    elseif(n>(p-1))
        dB(i,:) = 0.5*p*(B(i-1,:));
    else
        dB(i,:) = 0.5*p*(B(i-1,:)-B(i,:));
    end
end