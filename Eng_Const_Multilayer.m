%% ENGINEERING CONSTANTS OF A MULTILAYERED COMPOSITE LAMINATE ACCORDING TO CLT
% Author: Augusto Lima
% Politècnico di Milano, April 2024
clc
clear all
%% CLT assumptions
% 1) Each layer is quasihomogeneous and orthotropic
% 2) The laminate is thin with its lateral dimensions much larger than its
% thickness + in-plane loading -> plane stress state
% 3) All displacements are small wrt the thickness of the laminate
% 4) Displacements are continuous throught the laminate
% 5) In-plane displacements in the x and y-directions are linear fcns of z
% 6) Trasnverse shear are negligible (straight sections remain straight
% after deformation)
% 7) Elastic material behavior
% 8) Transverse normal strain is negligible

%% INPUTS
% ENGINEERING CONSTANTS - LAYER REFERENCE SYSTEM
E1 = 146000; %[GPa]
E2 = 8220; %[GPa]
nu12 = 0.32;
G12 = 4500; %[GPa]
S = [  1/E1,     -nu12/E1,  0;
     -nu12/E2,    1/E2,    0;
        0,           0,   1/G12]; %Compliance matrix

% ORIENTATION, NUMBER OF LAYERS, SEQUENCE OF LAYUP AND THICKNESSES
n = 8; % Number of layers
Orientation =[ 0; %First layer Orientation [deg]
               +45; %Second layer Orientation [deg] ...
               -45;
               90;
               90;
               -45;
               45;
               0]*pi/180 ;
Thicknesses =[ 0.15; %First layer Thickness [mm]
               0.15;  %Second layer Thickness [mm]
               0.15; 
               0.15;
               0.15;
               0.15;
               0.15;
               0.15] ; 
%% BUILDING LAYUP
% Converting all compliance/stiffness tensors of each layer to global
% reference system
Sxy = cell(n,1);
Qxy = cell(n,1);
for i = 1:n
    theta = Orientation(i);
    cosine = cos(theta);
    sine = sin(theta);
    T = [cosine^2  sine^2  2*cosine*sine;
     sine^2  cosine^2 -2*cosine*sine;
    -cosine*sine  cosine*sine  cosine^2-sine^2]; %Rotation matrix
    Sxy{i} = T'*S*T; %Compliance matrix
    Qxy{i} = inv(Sxy{i}); %Stiffness matrix
end

% Defining Midplane
t = sum(Thicknesses, "all");
z0 = t/2;

% Defining z-distances to local midplanes
t_rest = 0;
z = zeros(n,1);
for i = 1:n
    z_bar = t-t_rest-(Thicknesses(i)/2); %Starting from the bottom, this will give the location of the midplane
    z(i) = z_bar-z0; %This will position the z wrt the global midplane
    t_rest = t_rest + Thicknesses(i);
end

% Defining Layup cells
Layup = {Orientation, Thicknesses, z, Sxy, Qxy};


%% BUILDING COMPOSITE MATRICES
% Building matrices ABD
% A matrix
A = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:n
            Qk = Qxy{k};
            hk = z(k)+Thicknesses(k)/2;
            hk_  = z(k)-Thicknesses(k)/2;
            A(i,j) = A(i,j) + Qk(i,j)*(hk-hk_);
        end
    end
end
%B matrix
B = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:n
            Qk = Qxy{k};
            hk = z(k)+Thicknesses(k)/2;
            hk_  = z(k)-Thicknesses(k)/2;
            B(i,j) = B(i,j) + Qk(i,j)*(hk^2-hk_^2)/2;
        end
        if B(i,j) <= 10^(-12) %Just a minor adjustment
            B(i,j) =0;
        end
    end
end
%D matrix
D = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:n
            Qk = Qxy{k};
            hk = z(k)+Thicknesses(k)/2;
            hk_  = z(k)-Thicknesses(k)/2;
            D(i,j) = D(i,j) + Qk(i,j)*(hk^3-hk_^3)/3;
        end
    end
end
Q_composite = [A,B;B,D];

% Building matrices abc
A_ = inv(A);
B_star = -A_*B;
C_star = B*A_;
D_star = D-(B*A_)*B;

d = inv(D_star);
a = A_-(B_star*d)*C_star;
b = B_star*d;
c = -d*C_star;

q_composite = [a,b;c,d];


%% OBTAINING THE ENGINEERING CONSTANTS
Ex = 1/(t*a(1,1));
Ey = 1/(t*a(2,2));
Gxy = 1/(t*a(3,3));
nuxy = -a(2,1)/a(1,1);
nuyx = -a(1,2)/a(2,2);
etasx = a(1,3)/a(3,3);
etaxs = a(3,1)/a(1,1);
nuys = a(3,2)/a(2,2);
nusy = a(2,3)/a(3,3);









