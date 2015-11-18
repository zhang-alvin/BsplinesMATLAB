function [approx,x_coord]=bsplSetupAndSolve(p,N,npts,fun)
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

figure(100)
plot(x_coord,bspl);

%Get quadrature points
intorder = 5
[q,w] = lgwt(intorder,-1,1);
q = flipud(q);
w = flipud(w);

%Define the shape functions once
Nshp = zeros(nshp_l,length(q),eID);

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
            temp = 0;
            for k=1:length(q)
                temp=temp+Nshp(i,k,eID)*Nshp(j,k,eID)*w(k);
            end
            temp = temp*J(eID);
            K(ien(eID,i),ien(eID,j)) = K(ien(eID,i),ien(eID,j)) + temp;
        end
    end  
end

%RHS

F = zeros(nshp_g,1);
for eID = 1:N
    for i=1:nshp_l
        temp=0;
        for k=1:length(q)
            x = x_mesh(eID)*(1-q(k))/2+x_mesh(eID+1)*(1+q(k))/2;
            temp = temp+Nshp(i,k,eID)*fun(x)*w(k);
        end
        temp = temp*J(eID);
        F(ien(eID,i))=F(ien(eID,i))+temp;
    end
end

%Solve
coef = K\F;

%Solution and post-processing
approx=coef'*bspl;