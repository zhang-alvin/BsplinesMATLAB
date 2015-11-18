%% 11.14.15
%% Test to see how B-splines behave in different scenarios
% clear
% close all
% clc

%Number of plotting points
npts = 100;

%Desired function to project
%fun = @(x) x.^2;
%fun = @(x) exp(x);
%fun = @(x) sin(x*pi);
%fun = @(x) cos(x*2*pi);
%fun = @(x) tan(x);
alpha = 20;
fun = @(x) (1-exp(-alpha*x))/(1-exp(-alpha));

% numel_case = 4;
% nump_case = 2;
% err_tab = zeros(nump_case,numel_case);
% comp_tab = zeros(nump_case,numel_case);
% p = [3,10];
% for i=1:numel_case
%     N=2^(i-1);
%     for j=1:nump_case
%         [approx,x_coord,K,F,coef,L2,comp_time]=bsplSetupAndSolve(p(j),N,npts,fun,tic);
%         err_tab(j,i) = L2;
%         comp_tab(j,i) = comp_time;
%     end
% end

N = 32;
p = 3;

N = 1;
p = 10;

[approx,x_coord,K,F,coef,L2,comp_time]=bsplSetupAndSolve(p,N,npts,fun,tic);
figure(1)
plot(x_coord,fun(x_coord),'k',x_coord,approx,'o')