%% 11.14.15
%% Test to see how B-splines behave in different scenarios
clear
close all
clc

% Single element case



N = 2;
p = 2;
npts = 100;
fun = @(x) exp(x);
%fun = @(x) sin(x*pi);
%fun = @(x) cos(x*2*pi);
[approx,x_coord]=bsplSetupAndSolve(p,N,npts,fun);
figure(1)
plot(x_coord,approx,x_coord,fun(x_coord))