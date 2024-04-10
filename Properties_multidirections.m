% Script to calculate UD lamina composite's (under in-plane loading) 
% material properties for various loading directions

close all
clear
clc

%% Assumptions
% In-Plane Loading
% Orthotropic material
% Theta is defined as being positive in the counterclockwise direction
% starting from the x (loading) axis to principal axis 1

%% INPUTS

E1 =15.8; % [GPa]
E2 =12.9; % [GPa]
G12 =2.7; % [GPa]
nu12 =0.16; % 
theta_min =0; % [deg]
theta_max =90; % [deg]
theta_step = 1; % [deg]

%% Computing compliance matrix wrt principal axes

S = [1/E1     -nu12/E1   0;
     -nu12/E1 1/E2       0;
       0        0      1/G12];

properties = [];

%% Begin iteration
for theta = theta_min:theta_step:theta_max

% Computing transformation matrix
theta_rad = theta*pi/180;
m = cos(theta_rad);
n = sin(theta_rad);

T = [m^2  n^2  2*m*n;
     n^2  m^2 -2*m*n;
    -m*n  m*n  m^2-n^2];

% Comuting compliance matrix wrt loading ref frame

Sxy = T'*S*T;

% Computing Material properties wrt loading ref frame

Ex = 1/Sxy(1,1); % [GPa]
Ey = 1/Sxy(2,2); % [GPa]
Gxy = 1/Sxy(3,3); % [GPa]

nuxy = -Ex*Sxy(2,1); % 
etaxs = Ex*Sxy(3,1); %
etays = Ey*Sxy(3,2); %
etasx = Gxy*Sxy(1,3); %
etasy = Gxy*Sxy(2,3); %

mat = [theta, Ex, Ey, Gxy, nuxy, etaxs, etays, etasx, etasy];

properties = [properties; mat];

end

theta = properties(:,1);
Ex = properties(:,2);
Ey = properties(:,3);
Gxy =  properties(:,4);
nuxy = properties(:,5); % 
etaxs = properties(:,6); %
etays = properties(:,7); %
etasx = properties(:,8); %
etasy = properties(:,9); %

plot(theta,etasx)
