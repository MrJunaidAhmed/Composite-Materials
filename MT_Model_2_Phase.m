clc;
clear;
close all;

%% 
%Spherical inclusions: Aspect ratio = 1
%Volume fractions of inclusion nu1: 10%, 20%, 30%, 40%, 50%
%Stiffness contrast E1/E0: 100, 50, 20, 10, 5, 0

%Matrix Properties
E_0 = 3.42e9;                       %Young's modulus of matrix [Pa]
nu_0 = 0.32;                        %Poisson's ratio of matrix

G_0 = E_0 / (2*(1+nu_0));
K_0 = E_0 / (3*(1-2*nu_0));

%Eshelby's tensor volumetric and deviatoric parts:
K_S = (1 + nu_0) / (9*(1 - nu_0));
G_S = (4 - 5*nu_0) / (15*(1 - nu_0));

%Reinforcement Properties
nu_1 = 0.21;                        %Poisson's ratio of inclusion
v_1 = linspace(0, 0.5, 6);          %inclusion volume fraction
SC = [100, 50, 25, 10, 5, 0];       %stiffness constrast

A = length(SC);
B = length(v_1);

E_1 = zeros(A , 1);
G_1 = zeros(A , 1);
K_1 = zeros(A , 1);
G_hom_ratio = zeros(A , B);
K_hom_ratio = zeros(A , B);

for i = 1 : A
    E_1(i) = SC(i)*E_0;
    G_1(i) = E_1(i) / (2*(1+nu_1));
    K_1(i) = E_1(i) / (3*(1-2*nu_1));
    
    for j = 1 : B        
        %homogenized values normalized w.r.t. matrix values
        K_hom_ratio(j,i) = ((v_1(j)*K_1(i))/K_0 + (1 - v_1(j))*(1 + 3*K_S*((K_1(i)/K_0) - 1) )) / (v_1(j) + (1 - v_1(j))*(1 + 3*K_S*((K_1(i)/K_0) - 1) ));
        G_hom_ratio(j,i) = ((v_1(j)*G_1(i))/G_0 + (1 - v_1(j))*(1 + 2*G_S*((G_1(i)/G_0) - 1) )) / (v_1(j) + (1 - v_1(j))*(1 + 2*G_S*((G_1(i)/G_0) - 1) ));
    end
end

%Plotting the results

%figure(1)
%plot(v_1*100,G_hom_ratio)
%xlabel('Inclusions Volume Fraction %')
%ylabel('Normalized Effective Shear Modulus')
%legend({'E1/E0=100','E1/E0=50','E1/E0=20','E1/E0=10','E1/E0=5','E1/E0=0'})
%grid on

%figure(2)
%plot(v_1*100,K_hom_ratio)
%xlabel('Inclusions Volume Fraction %')
%ylabel('Normalized Effective Bulk Modulus')
%legend({'E1/E0=100','E1/E0=50','E1/E0=20','E1/E0=10','E1/E0=5','E1/E0=0'})
%grid on

%K_hom_ratio
%G_hom_ratio

A = [4, 3, 2, 1];
A(1, 1) = 7;
A

v_1 = linspace(0.1, 0.5, 5)

