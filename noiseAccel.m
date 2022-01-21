clc
%% Constants
sigma = 5.67e-8;                %stephan boltzmann constant (W*m^-2*K^-4)
c  = 3e8;                       %speed of light (m/s)
k_B = 1.38064852e-23;           %Boltzmann constant (m2 kg s-2 K-1)
%% Inputs
%Frequency Band
Nd = 1;
Nf = 1000;
f = logspace(-4, 0, Nf)';
f = f(:, ones(1, Nd));
vec = ones(1,Nf);   %vector
%How big can the dT (K) be in order to be within the noise budget?
%ltpda toolbox




%TM Characteristics
d   = 0.001;                    %gap (m)
s   = .030;                     %TM side length (m)
A   = s^2;                      %Area of face segment (1-D) (m^2)
A_wall = (s+2*d)^2;
rho = 2.0e4;                    %AuPt density (kg/m^3)
m   = rho*A*s;                  %TM mass (kg)
I   = 1/6*m*s^2;                %Moment of Inertia

%LISA TM Characteristics
d_LISA = 0.0032;
s_LISA = 0.046;
A_LISA = s_LISA^2;
rho_LISA = 2.0e4;
m_LISA = rho_LISA*A_LISA*s_LISA;
I = 1/6*m_LISA*s_LISA^2;

massrat = m_LISA/m;




%Thermal Characteristics
theta = 3e4;                    %activation temperature (K), in range [3e3, 3e4] (worst case 3e3)
Io = 0.5e-6;                     %kg/s^3, outgassing rate taken from accelNoise_geo_R0.m.         

%ask about ranging of the ADC from the Monday meeting
dT = 5e-3;
dTf = dT.*(1e-3./f).^(1/3);

P = 1e-5;                       %Reservoir Pressure (Pa)
T0 = 293;                       %Reservoir Temp (K)

% Brownian Residual Gas effects
M_gas = 18.01528;
m_gas = 2.9915076e-26;          %Reservoir Gas avg. mass (gas is water) (kg)

%% On-Board Magnetic Field Effects
X0 = 17e-6; mu0 = 1.25664e-6; sigma0 = 3.3e6; B = 1e-4; M = 1; gradBx = 1e-12;
Xtm = X0 - j*s^2*mu0*sigma0*pi.*f/12;

% a_obmf = (Xtm.*s^3/mu0*B+M)*gradBx;  %on board magnetic field
% a_mff = 3^(-1/2)*abs(Xtm)*s^3/(m*mu0)*deltaB*partialB; %magnetic field fluctuation
% a_mfgf = 3^(-1/2)/m*(s^3*X0*B/mu0 + mu_x)*partialdeltaB; %magnetic field gradient fluctuation
% a_ACmf = 5*s^2*Xtm/(m*mu0)*BmaxAC*partialBAC; %Down-Converted AC Magnetic Field


%Unknowns: term name (description [value for LISA, absent if not listed])
%        : X0 (magnetic susceptibility at f=0 Hz
%        : sigma0 (electrical conductance of proof mass [3.3×106 A/V/m]) material property?
%        : B (DC magnetic field (maximum component) [5 ?T])
%        : gradBx (magnetic field difference on each side of test mass [1.2566×10-6 N/A2])
%        : M (remnant magnetisation)
%        : deltaB (DC magnetic field gradient (maximum component) [1 ?T/m])
%        : partialB (magnetic field fluctuation [T/?Hz])
%        : mu_x (proof mass permanent magnetic moment [20×10-9 Am2])
%        : partialdeltaB (magnetic field gradient fluctuation [T/m/?Hz])
%        : BmaxAC (AC magnetic field maximum [0.5 ?T])
%        : partialBAC (Magnetic field fluctuation above MBW [T/?Hz])

%% Earth Magnetic Field effects

% A lot of unknowns here. Use SPENVIS to get the magnetic field
% fluctuations.
% Includes Lorentz force and non-lorentz forces
%% Electrical Forces

%Electrical Force on Proof Mass from Fluctuating Charge

%%%%%%%%%%% Stray DC Voltages and TM Charge Effects %%%%%%%%%%%%
% 1. Weber, Possibilities for measurement and compensation of stray DC electric fields acting on drag-free test masses
% 2. Charge-induced force-noise on free-falling test masses: results from LISA Pathfinder
% 3. Pollack, Temporal Extent of Surface Potentials between Closely Spaced Metals
% 4. Weber, Characterization of disturbance sources for LISA: torsion pendulum results
deltax = 200e-3;                        % potential of electrodes
lamda_eff = 260;                        % e/s, Weber

delV = 100e-3;
del = 10^-5;
q_e = 1.602e-19; %C, Elementary Charge

%a_Qm = sqrt(lamda)*q_e*partV_eQ/(sqrt(2)*m*pi*f*l_0); %Electrical Force on Proof Mass from Fluctuating Charge
%a_partV = 2*Cx/m/d*sqrt(partV_e^2+(Q_TM/C_g)^2)*partV; %Electrical Force on Proof Mass from Fluctuating Voltage
%a_Nyquist = 2*Q_TM/(m*d*(2*C_x+C_g))*sqrt(k_B*T*C_x*part_C/(pi.*f)); %Nyquist Noise Converted into Brownian Noise
%Electrostatic Readout Back-Action had an effect of 5e-18 and was
%considered negligible
%Electrical Force on Proof Mass due to Digital Resolution of Actuation
%Signal also extremely small
%a_T_DC_phi = 1/m*sqrt(2*C_x0/d*T_phi_DC/r_x)*partU_act100; %ctrical Force Noise from DC Torque actuation in ?
%%%Unknowns%%%

%For a_Qm:
%lamda charge rate in number of elementary charges per s [500 s-1]
%partV_eQ stray DC electrode potential for charge coupling [10 mV]
%l_0 charging length [44 mm]

%For a_partV:
%C_x capacitance on sensitive axis [2.31 pF]
%Q_TM proof mass charge [1.6 pC]
%?Ve stray DC voltage potential [0.1 V]
%Cg total capacitance of proof mass to ground [20.99 pF]
%?V Quasi-DC voltage fluctuations at electrodes [1?10-6 V/?Hz]

%For a_Nyquist:
%?C sensing capacity loss angle [1 ?rad]

%For a_T_DC_phi:
%?Uact,100Hz actuation voltage fluctuation at 100 Hz [10 ?V/?Hz]
%T?,DC maximum DC torque around ? in science mode [7?10-13 Nm]
%Cx0 actuation electrode capacitance in x [2.31 pF]
%rx actuation electrode lever arm in x [11.15 mm]
%TGRS electrode housing temperature [293 K]
%kB Boltzmann constant [1.38062×10-23 J/K]
%?sens sensing capacity loss angle [1 ?rad]


a_charge = 7e-15*(deltax/200e-3)*(lamda_eff/260)^(1/2)*(0.1e-3./f)*massrat;
%dielectric noise
a_diel = 5e-16*(delV/100e-3)*(del/10^-5)^(1/2)*(0.1./f).^(1/2)*massrat;

figure
loglog(f,a_charge,f,a_diel)
grid()
legend('Charge','Dielectric')
title('Thermal effects')
xlabel('Frequency (Hz)')
ylabel('Acceleration noise (m/s^2/Hz^{1/2})')
%where we have defined the translational DC bias

%% Brownian Noise

%%%%%%%%%%% Residual Gas Damping Effects %%%%%%%%%%%%%%5           
% 1. Cavalleri, PhysRevLett.103.140601
% 2. Cavalleri, Gas damping force noise on a macroscopic test body in an infinite gas reservoir
B_inf_tr = P*s^2*(1+pi/8)*(32*m_gas/pi/k_B/T0)^(1/2);           %Brownian damping coefficient for cube translating
B_tr = B_inf_tr/log(s/d)/(d/s)^2;

S_F = 4*k_B*T0*real(B_tr);

%a_diel = 2*C_sens/m/d*(partV_e+Qmax/Ctot)*sqrt(4*T0*kB*part_sens/(2*pi.*f*C_sens)); %Dielectric Losses
%a_gas = 1e-15*sqrt(P/1e-6); %Residual Gas Impacts on Proof Mass
%a_eddy = (3/2^0.5*pi)^(1/3)*sqrt(kB*T0*s^5*sigma0/5)*gradB/m
%a_ferro = 2*max(gradB)/pi/m*sqrt(kB*T0*s^3*partX/mu0./f); %Magnetic
%Impurities / Magnetic Viscosity in the Ferromagnetic Component, 
%ASD page 51 for more info
a_rgi = sqrt(S_F)/m*ones(size(f(:, 1))); %same as a_gas, check if same values. Likely won't scale accordingly because of geometry difference

%%%Unknowns%%%
%?Ve stray DC electrode potential [100 mV]
%Csens Single electrode sensing capacity [1.15 pF]
%Qmax proof mass charge [1.6 pC]
%?sens sensing capacity loss angle [1 ?rad]
%?0 proof mass electrical conductance [3.3 MA/V/m]
%?B DC magnetic field gradient (maximum element) [1 ?T/m]

%% Thermal Effects
 
%%%%%%%%%%% Conductance Calculations %%%%%%%%%%%%%
%These calculations come from: https://uspas.fnal.gov/materials/15ODU/Session1_Fundamentals.pdf
%These calculations for gas conductance are necessary to come up with a
%ratio (C_L2G) to relate the LISA conductance channels to the GRACE
%channels. Each variable is given a suffix of "g" for GRACE and "l" for 
%LISA.
L_g = (s+d)/d;
L_l = (s_LISA+d_LISA)/d_LISA;
R_g = s/d;
R_l = s_LISA/d_LISA;
K_g = 1.0663+0.0471*R_g-0.000848*R_g^2;
K_l = 1.0663+0.0471*R_l-0.000848*R_l^2;
D_g = 0.7749+0.1163*R_g-0.000372*R_g^2;
D_l = 0.7749+0.1163*R_l-0.000372*R_l^2;
E_g = 0.6837+0.02705*R_g-0.000616*R_g^2;
E_l = 0.6837+0.02705*R_l-0.000616*R_l^2;
a_s_g = 1/(1+L_g*(1+R_g)/(2*R_g));
a_s_l = 1/(1+L_l*(1+R_l)/(2*R_l));
a_l_g = 1/(1+3*L_g*(1+R_g)/(8*R_g*K_g));
a_l_l = 1/(1+3*L_l*(1+R_l)/(8*R_l*K_l));
a_g = (1+D_g*(L_g/R_g)^E_g)/(1/a_s_g+(D_g*(L_g/R_g)^E_g)/a_l_g);
a_l = (1+D_l*(L_l/R_l)^E_l)/(1/a_s_l+(D_l*(L_l/R_l)^E_l)/a_l_l);
C_g = a_g*3.64*d*(s+d)*(T0/M_gas)^(1/2);
C_l = a_l*3.64*d_LISA*(s_LISA+d_LISA)*(T0/M_gas)^(1/2);
C_L2G = C_g/C_l;        %Conductance transfer functio from Lisa to Grace

%Conductance of the hole
d_hole = 0.0035;
l_hole = 0.01;
v_hole = sqrt(8*k_B*T0/(pi*m_gas));
a_hole = 1.0667/(1.0667+(l_hole/d_hole)^0.94629); %transmission probability fit for l/d < 50
Corif = a_hole*v_hole/4*d_hole^2/4*pi;

% 1. Carbone PhysRevD.76.102003
% 2. R Dolesi et al 2003 Class. Quantum Grav. 20 S99
% 3. Astruim LISA Requirement Breakdown, LISA-ASD-TN-5001, 29 Nov 2009

Chole = 4.3e-3;             %m^3/s, vent hole flow, taken from accelNoise_geo_R0.m.
Cch = C_g;
%Cch = (1.33+1.4)*Chole;     %Channel Conductance ratio from LISA,  ASD
%Cch = Cch*C_L2G;            %Scaling between LISA and GRACE is applied
Ceff = 2*Cch+Chole;         %Effective Conductance,  R Dolesi et al 2003 Class. Quantum Grav. 20 S99

dF_RdT = A*P/4/T0; %N/K
dF_RPdT = A*8/3*sigma/c*T0^3; %N/K
dF_OGdT = 5.8e-13; %A_wall*A*Io*theta/(T0^2)/Ceff; %N/K


a_R = dF_RdT.*dTf/m;
a_RP = dF_RPdT.*dTf/m;
a_OG = dF_OGdT.*dTf/m;
dFdT = dF_RdT+dF_RPdT+dF_OGdT;
uppLim = dFdT*5e-3/1.96

a_TG = sqrt(a_R.^2 + a_RP.^2 + a_OG.^2);


figure
loglog(f,a_TG,f,a_R,f,a_RP,f,a_OG)
grid()
legend('Total Thermal Effects','Radiometer','Radiation Pressure','Outgassing')
xlabel('Frequency (Hz)')
ylabel('Acceleration Noise (m/s^2/Hz^{1/2})')
%% Spacecraft Self Gravity

%% Others
Chole = 4.3e-3;             %m^3/s, vent hole flow, taken from accelNoise_geo_R0.m.
Cch = C_g;
%Cch = (1.33+1.4)*Chole;     %Channel Conductance ratio from LISA,  ASD
%Cch = Cch*C_L2G;            %Scaling between LISA and GRACE is applied
Ceff = 2*Cch+Chole;         %Effective Conductance,  R Dolesi et al 2003 Class. Quantum Grav. 20 S99

a_R_lpf = 0.046^2*1e-6/(2*293*1.96).*dTf;
a_RP_lpf = 8*1*sigma*293^3*0.046^2/(3*1.96*3e8).*dTf;
a_OG_lpf_1 = 2.5115e-12.*dTf/1.96; %0.046^2*0.0028*0.5e-6*30000/(1.96*293^2*0.0043*(1+2*(1.33+1.4))).*dTf;
a_OG_lpf_2 = 0.046^2*0.0028*0.5e-6*30000/(1.96*293^2*0.0043*(1+2*(1.33+1.4))).*dTf;


a_TG = a_R_lpf + a_RP_lpf + a_OG_lpf_2;

figure
loglog(f,a_R_lpf,f,a_RP_lpf,f,a_OG_lpf_1)
grid()
legend('Radiometer','Radiation Pressure','Outgassing')
title('Thermal effects LPF Molflow Simulation')
xlabel('Frequency (Hz)')
ylabel('Acceleration noise (m/s^2/Hz^{1/2})')

figure
loglog(f,a_R_lpf,f,a_RP_lpf,f,a_OG_lpf_2)
grid()
legend('Radiometer','Radiation Pressure','Outgassing')
title('Thermal effects LPF Equation')
xlabel('Frequency (Hz)')
ylabel('Acceleration noise (m/s^2/Hz^{1/2})')


%% Stray DC voltages and TM charge

% charge noise*stray field &&& noise in stray field*charge 2 largest ones probably
% what's avg tM charge, noise in charge (comes from discharging system),
% stray electric field in the sensor, and offsets not perfecty stable. Well
% understood for LISA

% Quick rough estimate for paper is good, more detailed later maybe? 
% Hatch potentials, force disturbance
% In LISA, gaps are big enough that this is not big effect
% Find an expression for forces and noise in terms of capacitances and gaps
% (capacitance gradients)

% In sumner paper: Equation 22 is the biggest contributer. Q/CT is a
% voltage, can probably get from LISA data. Use same assumptions as we did
% in LISA, 70 mV
% Svi we can use LISA number (noise in voltage) 200 microV/Hz^1/2
% Ci - capacitance of x 
% gradient is just derivative of capacitance, E0*A/d^2 (divided by d
% again?)


% Sa8 we can calculate. Sq is the noise that comes from discharging TM
% coontinuously
% Vi is the stray potential. 200 microV/Hz^1/2
% 8&11 almost same equation, noise of Q in 8 and noise of Vi in 11

% Vi might not need DC compensation, could be 100 instead of 3

%SQ^1/2 lamda_eff is 5000 to calculate charge noise
E0 = 8.854e-12;
A = (32e-3)^2; %mm^2, area of TM face
d = 1e-3;
Ci = E0*A/d;
dCdK = Ci/d;
C_T = 6*Ci;

e = 1.602e-19; % Coulombs, charge of an electron 
M = 1.96; % kg, mass of proof mass (TM)
lamda_eff = 5000; % e/s
Vtm = 70e-3; %70 mV, Q/C_T

%TM Charge
S_Q = e./(2*pi.*f).*sqrt(2*lamda_eff);

%Stray Potentials
Vi = 0.1; % 100 mV
S_Vi = sqrt(8.83e-16./(f.^2) +  2.73e-13./f)*Vi;


Sa8 = 1/(M*C_T)*Vi*dCdK*S_Q;
Sa11 = Vtm/M*dCdK*S_Vi;

figure
loglog(f,Sa11,f,Sa8,f,sqrt(Sa11.^2+Sa8.^2))
xlabel('Frequency (Hz)')
ylabel('Acceleration Noise (m/s^2/Hz^1/2)')
grid()

%% Total Acceleration Noise
a_th = a_R+a_RP+a_OG;
% a_elec = a_charge+a_diel;
a_charge = Sa8+Sa11;
a = a_th+a_rgi+a_charge;
a_budget = 3e-12.*sqrt(1+(f./.008).^4);

figure
% loglog(f,a,f,a_th,f,a_rgi,f,a_elec,f,a_budget,f1,sqrt(pxx1))
loglog(f,a,f,a_th,f,a_charge)
legend('Total Noise','Thermal Noise','Charge related noise')
xlabel('Frequency (Hz)')
ylabel('Acceleration Noise (m/s^2/Hz^{1/2})')
grid()

