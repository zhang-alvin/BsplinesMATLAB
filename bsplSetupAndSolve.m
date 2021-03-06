function [approx,x_coord,K,F,coef,L2,comp_time]=bsplSetupAndSolve(p,N,npts,fun,before)
%% function that sets up BFEM and solves for the resulting solution

% Single element case

cont = p-1;
mult = p-cont;
nshp_g = p+1+(N-1)*(mult);
nshp_l = p+1;
x_mesh = linspace(0,1,N+1);

%Get Jacobian
J = zeros(N,1);
for i=1:N
    J(i) = x_mesh(i+1)-x_mesh(i);
end

%Get knot vector in parametric space
numknots = (N-1)*(p-cont)+2+2*p;
knot = zeros(numknots,1);
knot(1:p+1) = 1;
knot(length(knot)-p:length(knot)) = length(x_mesh);
count = 2;
for i=p+2:length(knot)-p-1
    knot(i) = count;
    count = count+1;
end

%Get Extraction Operator
if(p==1)
    C = zeros(nshp_l,nshp_l,N);
    for i=1:N
        C(:,:,i) = eye(nshp_l,nshp_l);
    end
else
    C = Bezier_extract(knot);
end

%Get IEN Array with only p-1 splines

%ien[eID][shpfun#]

ien = zeros(N,nshp_l);
count = 1; %following matlab convention
for i=1:N
    for j=1:nshp_l
        ien(i,j) = count;
        count = count+1;
    end
    count = count - p;
end

%Get Splines
x_spline = linspace(-1,1,npts);
bspl = zeros(nshp_g,npts*N);
for i = 1:nshp_l
    for eID = 1:N
        index = ien(eID,i);
        bspl(index,1+npts*(eID-1):npts*(eID)) = getBernstein(p,x_spline,'spline',C(i,:,eID)');
    end
end

%Get plotting coordinates
x_coord = zeros(npts*N,1);
for eID=1:N
    x_coord(1+npts*(eID-1):npts*(eID)) = x_mesh(eID)*(1-x_spline)/2 + x_mesh(eID+1)*(x_spline+1)/2;
end

%figure(100)
%plot(x_coord,bspl);

%Get quadrature points
intorder = 2*(p+1);
[q,w] = lgwt(intorder,-1,1);
q = flipud(q);
w = flipud(w);

%Define the shape functions once
Nshp = zeros(nshp_l,length(q),N);

for eID=1:N
    for i=1:nshp_l
        Nshp(i,:,eID) = getBernstein(p,q','spline',C(i,:,eID)');
    end
end

%LHS

K = zeros(nshp_g);

for eID = 1:N
    for i=1:nshp_l
        for j=1:nshp_l
            temp = Nshp(i,:,eID).*Nshp(j,:,eID)*w;
            temp = temp*J(eID);
            K(ien(eID,i),ien(eID,j)) = K(ien(eID,i),ien(eID,j)) + temp;
        end
    end  
end

%RHS

F = zeros(nshp_g,1);
for eID = 1:N
    for i=1:nshp_l
        x = x_mesh(eID)*(1-q)/2+x_mesh(eID+1)*(1+q)/2;
        temp = Nshp(i,:,eID).*fun(x)'*w;
        temp = temp*J(eID);
        F(ien(eID,i))=F(ien(eID,i))+temp;
    end
end

%Apply Interpolatory BC
F(1) = fun(min(x_coord));
F(end) = fun(max(x_coord));
K(1,1) = 1;
K(end,end) = 1;
K(1,2:end)=0;
K(end,1:end-1)=0;
for i=2:nshp_g-1
    F(i) = F(i) - K(i,1)*F(1);
    K(i,1) = 0;
end
for i=nshp_g-1:-1:2
    F(i) = F(i) - K(i,end)*F(end);
    K(i,end)=0;
end

%Solve
coef = K\F;

comp_time = toc(before);

%Solution and post-processing
approx=coef'*bspl;

%Compute the L2 error

L2=0.0;
for eID = 1:N
    x = x_mesh(eID)*(1-q)/2+x_mesh(eID+1)*(1+q)/2;
    L2=L2+w'*(Nshp(:,:,eID)'*coef(ien(eID,:))-fun(x)).^2;
end
L2 = sqrt(L2);