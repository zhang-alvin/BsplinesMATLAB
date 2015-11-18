function [B]=getBernstein(p,xi,mode,beta)

%Default
if ~exist('mode','var')
    mode = 'bernstein';
end

numqpt = length(xi);
if strcmp(mode,'bernstein')
    B  = zeros(p+1,numqpt);
    for i=1:p+1
        n = i-1;
        B(i,:) = nchoosek(p,n)*(1-0.5*(xi+1)).^(p-n).*(0.5*(xi+1)).^n;
    end
elseif strcmp(mode,'spline')
    %B  = zeros(numqpt,1);
    off = -1;
    J = 0.5;
    B = DeCastelJau(p,beta,xi,J,off);
else
    error('no proper mode specified. use either bernstein or spline');
end

end %function end


