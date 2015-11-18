function [value] = DeCastelJau(p,beta,t0,J,off)
%% p is polynomial order
%% betas is the coefficients associated with each polynomial Bernstein basis
    % Note that beta is in column major form so that beta for the first 2nd
    % order bernstein polynomial is beta = [1;0;0] and not beta = [1,0,0];
%% t0 is the desired point
%% J is a scaling factor
%% off is the offset from 0

B0 = inline('1-J*(t-off)','t','J','off'); %1-t
B1 = inline('J*(t-off)','t','J','off'); %t

numpts = length(t0);

beta_test = zeros(p+1,p,numpts);
beta_test(:,1,:) = repmat(beta,1,numpts);  
if(numpts==1)
    for j=2:p+1
        for i=1:p
            beta(i,j) = beta(i,j-1)*B0(t0,J,off)+beta(i+1,j-1)*B1(t0,J,off);
            %beta_test(i,j,:) = reshape(beta_test(i,j-1,:),1,p+1).*B0(t0,J,off)+reshape(beta_test(i+1,j-1,:),1,p+1).*B1(t0,J,off);
        end
    end
    value = beta(1,p+1);
else
    for j=2:p+1
        for i=1:p
            %beta(i,j) = beta(i,j-1)*B0(t0,J,off)+beta(i+1,j-1)*B1(t0,J,off);
            beta_test(i,j,:) = reshape(beta_test(i,j-1,:),1,numpts).*B0(t0,J,off)+reshape(beta_test(i+1,j-1,:),1,numpts).*B1(t0,J,off);
        end
    end
    value = reshape(beta_test(1,p+1,:),1,numpts);
end
