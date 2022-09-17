clc;
clear;
close all;

%% 
%Spherical inclusions: Aspect ratio = 1
%Volume fractions of inclusions v1 = v2: 0%, 10%, 20%, 30%, 40%, 50%
%Stiffness contrast for solid inclusions: E1/E0: 100, 50, 20, 10, 5

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
v_1 = linspace(0, 0.25, 6);          %inclusion volume fraction
B = length(v_1);

SC = [100, 50, 25, 10, 5];       %stiffness constrast of inclusion 1
A = length(SC);

E_1 = zeros(A , 1);
G_1 = zeros(A , 1);
K_1 = zeros(A , 1);
GHR1 = zeros(A , B);
KHR1 = zeros(A , B);

for i = 1 : A
    E_1(i) = SC(i)*E_0;
    G_1(i) = E_1(i) / (2*(1+nu_1));
    K_1(i) = E_1(i) / (3*(1-2*nu_1));

    for j = 1 : B        
        %homogenized values normalized w.r.t. matrix values
        KHR1(i,j) = ((2*v_1(j)*K_1(i))/K_0 + (1 - 2*v_1(j))*(1 + 3*K_S*((K_1(i)/K_0) - 1) )) / (2*v_1(j) + (1 - 2*v_1(j))*(1 + 3*K_S*((K_1(i)/K_0) - 1) ));
        GHR1(i,j) = ((2*v_1(j)*G_1(i))/G_0 + (1 - 2*v_1(j))*(1 + 2*G_S*((G_1(i)/G_0) - 1) )) / (2*v_1(j) + (1 - 2*v_1(j))*(1 + 2*G_S*((G_1(i)/G_0) - 1) ));
    end
end

%cavities
nu_2 = 0;
SC2 = 0;
v_2 = v_1;
C = length (v_2);
E_2 = SC2*E_0;
G_2 = E_2 / (2*(1+nu_2));
K_2 = E_2 / (3*(1-2*nu_2));
GHR2 = zeros(1 , C);
KHR2 = zeros(1 , C);

for j = 1 : C        
    %homogenized values normalized w.r.t. matrix values
    KHR2(1, j) = ((2*v_2(j)*K_2)/K_0 + (1 - 2*v_2(j))*(1 + 3*K_S*((K_2/K_0) - 1) )) / (2*v_2(j) + (1 - 2*v_2(j))*(1 + 3*K_S*((K_2/K_0) - 1) ));
    GHR2(1, j) = ((2*v_2(j)*G_2)/G_0 + (1 - 2*v_2(j))*(1 + 2*G_S*((G_2/G_0) - 1) )) / (2*v_2(j) + (1 - 2*v_2(j))*(1 + 2*G_S*((G_2/G_0) - 1) ));
end

%Voigt model
GHR = zeros(A , B);
KHR = zeros(A , B);

for i = 1 : A
   for  j = 1 : B
       GHR(i,j) = 0.5*GHR1(i,j) + 0.5*GHR2(1,j);
       KHR(i,j) = 0.5*KHR1(i,j) + 0.5*KHR2(1,j);
   end
end

%Plotting the results

figure(1)
plot(2*v_1*100, transpose(GHR))
xlabel('Inclusions Volume Fraction %')
ylabel('Normalized Effective Shear Modulus')
legend({'E1/E0=100','E1/E0=50','E1/E0=20','E1/E0=10','E1/E0=5'})
grid on

figure(2)
plot(2*v_1*100, transpose(KHR))
xlabel('Inclusions Volume Fraction %')
ylabel('Normalized Effective Bulk Modulus')
legend({'E1/E0=100','E1/E0=50','E1/E0=20','E1/E0=10','E1/E0=5'})
grid on