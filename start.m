clc
clear all

%% Constants

e = 1.602177e-19;               % electron charge (C)
mu0 = 1.2566e-6;                % mag vacuum permeability (N/A^2)
sigma = 5.67032e-8;             % Stefan-Boltzmann const. (W/m^2/K^-4)
c = 3e8;                        % speed of light (m/s)
G = 6.7e-11;                    % gravitational constant (m^3/kg/s^2)
GM = 398600e9;                  % std. grav. param. (m^3/s^2)
eps0 = 8.85e-12;                % vacuum permittivity (F/m)
kB = 1.38e-23;                  % boltzmann const. (J/K)
re = 6371e3;                    % Earth's diameter (m)

%% Parameters

Nsides = 6;
w_in = 100000*2*pi; %is it this in rad/s or just 100 kHz?

%% Formulas provided from Giacomo Ciani Thesis Section 2.2
f_str = 0;              %Stray forces
m = 1;                  %Test Mass mass (kg)
k = 0;                  %Global spring constant
dx = 0;                 %TM-spacecraft relative displacement

a_n = f_str/m-k*dx/m;    %acceleration noise

%% TM voltage

% Polarization due to surrounding voltages
V_tmj = 0;
for j=1:Nsides
    V_tmj = V_tmj + C_surf(j)/C_tot*Vi
end

% Injection Bias
V_tminj = alpha*V_inj*sin(w_in*t)

% TM Charge
V_tmQ = Q/C_tot


F_x = 1/2*sum(?C_i/?x*(V_tm-V_i)^2)