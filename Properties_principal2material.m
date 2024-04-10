% Script to calculate UD lamina composite's (under in-plane loading) 
% material properties given a loading direction

close all
clear
clc

%% Assumptions
% In-Plane Loading
% Orthotropic material
% Theta is defined as being positive in the counterclockwise direction
% starting from the x (loading) axis to principal axis 1

%% INPUTS

E1 =146; % [GPa]
E2 =8.22; % [GPa]
G12 =4.5; % [GPa]
nu12 =0.32; % 
theta =45; % [deg]

%% Computing compliance matrix wrt principal axes

S = [1/E1     -nu12/E1   0;
     -nu12/E1 1/E2       0;
       0        0      1/G12];

%% Computing transformation matrix
theta = theta*pi/180;
m = cos(theta);
n = sin(theta);

T = [m^2  n^2  2*m*n;
     n^2  m^2 -2*m*n;
    -m*n  m*n  m^2-n^2];

%% Comuting compliance matrix wrt loading ref frame

Sxy = T'*S*T;

%% Computing Material properties wrt loading ref frame

Ex = 1/Sxy(1,1); % [GPa]
Ey = 1/Sxy(2,2); % [GPa]
Gxy = 1/Sxy(3,3); % [GPa]

nuxy = -Ex*Sxy(2,1); % 
etaxs = Ex*Sxy(3,1); %
etays = Ey*Sxy(3,2); %













