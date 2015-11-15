function [coord,final_splines] = BsplineGenerator_GalerkinProjection(knot,resolution)

%Galerkin Projection of B-splines using Extraction Operator (part 1/3) 

%This script will generate the B-splines required for the given problem using 
%Cox-de Boor recursion. This can be used to generate B-splines for any open
%knot vector

%clear 
%clc
%close all

%% Declare the knot vector

%knot = [0,0,0,0,1,1,1,1,2,2,2,2]; %test case 1
%knot = [0,0,0,0,0,0,1,2,2,2,2,2,2]; %test case 2
%knot = [-0.05, -0.05, -0.05, -0.045, -0.036696, -0.0229047, 1.02109e-12,
%1.02109e-12, 1.02109e-12]; %scalar
%knot = [-0.0125 -0.0125 -0.0125 -0.0121 -0.011499 -0.0105959 -0.00923893 -0.0072 -0.0072 -0.0072]; %Couette/BL Case/F.S.
%resolution = 1000
%knot = [-0.0125 -0.0125 -0.0125 0 0 0]; %1 layer
%knot = [-0.0125 -0.0125 -0.0125 -0.00625 0 0 0]; %2 layers
%knot = [-0.0125 -0.0125 -0.0125 -0.009375 -0.00625 -0.003125 0 0 0]; %4 layers
%knot = [-0.0125 -0.0125 -0.0125 -0.0109375 -0.009375 -0.0078125 -0.00625 -0.0046875 -0.003125 -0.0015625 0 0 0];%8layers

[x,y] = mode(knot); %gets the second mode of the knot vector
p = y-1;
knotsize = length(knot);
num_poly = knotsize-(p+1);

%This defines the domain for the splines
x = linspace(knot(1),knot(knotsize),resolution);

%% Find indices along the domain vector,x, that correspond to each knot 

position_index = zeros(knotsize,1);
for i = 1:knotsize
    if i<=p+1 
        position_index(i) = find(x==knot(i));
    else
        position_index(i) = find(x>=knot(i),1);
    end
end

N = zeros(num_poly,p+1,length(x)); %initialize the Bspline matrix


%% p = 0; Define the base case - this will yield rectangular pulses in certain knot intervals

N_0_condition = knot(1:num_poly)<knot(2:num_poly+1); %If the knot interval is greater than 0, there will be a rectangular pulse
for i = 1:num_poly
    if N_0_condition(i) == 1
        N(i,1,position_index(i):(position_index(i+1)-1)) = 1; %The -1 is necessary because you only want to make the rectangular pulses within a single element
    end
    if i == num_poly %This is needed because otherwise the very last point will be 0 instead of 1.
        N(i,1,length(x))= 1;
    end
end

%% p = 1; Cox-de Boor recursion starts

 ksi = ones(1,1,length(x));
 ksi(1,1,:) = x; %This is needed because the matrix N has funny dimensions
 
 for k = 2:p+1 %k being 1 higher than it should be due to MATLAB indexing; this represents the polynomial order
     for i = 1:num_poly
         if i==num_poly %The num_polyth spline will always refer to a num_poly+1 spline that is undefined so it must be accounted for
             N(i,k,:)= (ksi-knot(i))./(knot(i+k-1)-knot(i)).*N(i,k-1,:); 
         elseif (knot(i+k-1)-knot(i) == 0) %If the denominator of a term ever goes to 0, that term must be ignored
             N(i,k,:) = (knot(i+k-1+1)-ksi)./(knot(i+k-1+1)-knot(i+1)).*N(i+1,k-1,:);
         elseif (knot(i+k-1+1)-knot(i+1) == 0)
             N(i,k,:)=(ksi-knot(i))./(knot(i+k-1)-knot(i)).*N(i,k-1,:);
         else
             N(i,k,:) = (ksi-knot(i))./(knot(i+k-1)-knot(i)).*N(i,k-1,:)+(knot(i+k-1+1)-ksi)./(knot(i+k-1+1)-knot(i+1)).*N(i+1,k-1,:);
         end
     end
 end
 
 %Store the final splines into an array for ease of plotting
 final_splines = zeros(num_poly,length(x));
 for i = 1:num_poly
    final_splines(i,:)=reshape(N(i,p+1,:),1,length(x));
 end
% line_types = ['-', ':', '', '.', 'd', '.k'];
%figure(1)
%  hold on
%  for i=1:num_poly
%  plot(x,final_splines(i,:),[line_types(i) 'k'],'LineWidth',2.5);
%plot(x,final_splines,'LineWidth',5);
%axis([knot(1) knot(knotsize) 0 1])
%  end
%  hold off
 
%legend('N_1_,_2', 'N_2_,_2', 'N_3_,_2','N_4_,_2','location','southoutside')

%title('Bsplines (p=1; knot=[0,0,1,2,2])')
 
% plot(x,final_splines(1,:),'k',x,final_splines(2,:),':k',x,final_splines(3,:),'.k',x,final_splines(4,:),'--k',x,final_splines(5,:),'dk',x,final_splines(6,:),'-.k','LineWidth',2.5);
% axis([knot(1) knot(knotsize) 0 1]);
% legend('N_1_,_2', 'N_2_,_2', 'N_3_,_2','N_4_,_2','N_5_,_2','N_6_,_2','orientation','horizontal','location','northoutside')
%  %title('B-splines for Boundary Layer Mesh Knot Vector')
% figure(2)
%GalerkinProjector_BezierExtraction
coord=x;