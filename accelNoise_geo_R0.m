% accelNoise_geo_R0.m
% j.w.conklin
% 2019.08.26
%
% references
% [D] Domenico Gerardi et al 2013
% [A] Astruim LISA Requirement Breakdown, LISA-ASD-TN-5001, 29 Nov 2009
% [S] Schumaker, CQG 2003
% [W] Weber et al, "Gas Damping ...", arXiv and PRL 2009
% [W2] Weber PRL 2012
% [R] R. Dolesi, LISA 9 
% [SMAD] Space mission analysis and design
% [NASA] NASA LISA requirements document

close all
clear all

% constants
e = 1.602177e-19;               % electron charge (C)
mu0 = 1.2566e-6;                % mag vacuum permeability (N/A^2)
sigma = 5.67032e-8;             % Stefan-Boltzmann const. (W/m^2/K^-4)
c = 3e8;                        % speed of light (m/s)
G = 6.7e-11;                    % gravitational constant (m^3/kg/s^2)
GM = 398600e9;                  % std. grav. param. (m^3/s^2)
eps0 = 8.85e-12;                % vacuum permittivity (F/m)
kB = 1.38e-23;                  % boltzmann const. (J/K)
re = 6371e3;                    % Earth's diameter (m)

% frequency band
Nd = 1;
Nf = 1000;
f = logspace(-4, 0, Nf)';
f = f(:, ones(1, Nd));

%%%%%%%%%%%%%%%   BEGIN PARAMETER INPUTS   %%%%%%%%%%%%%%%

% test mass design parameters
s = 0.030;                      % length (m)
d = 0.001;                      % gap (m)
rho = 0.9e4;                    % density (kg/m^3) [BeCu]
% rho = 2.0e4;                    % density (kg/m^3) [AuPt]
dd = 1e-6;                      % variations in gap size (m)
chi0 = 2e-6;                    % susceptibility at f = 0 [BeCu]
% chi0 = 2e-5;                    % susceptibility at f = 0 [AuPt]
dchi = 3e-7;                    % susceptibility, imag. comp.
Mr = 1e-8;                      % TM remnant mag. moment (Am^2) [R,PTB]
% Mr = 4e-8;                    % TM remnant mag. moment (Am^2) [A]
% Mr = 2e-8;                    % TM remnant mag. moment (Am^2) [D]
q = 1e7*e;                      % test mass charge (C) [R]
    % note on charge: GP-B = 0.1 mV/day, C = 5e-9 F -> 5e-13 C/day
sigma0 = 0; %3e6;               % elec. conductance (A/V/m) [BeCu,check]
Io = 0.5e-6;                    % outgassing rate (kg.sec^3) [A]
Thetao = 30000;                 % activation temperature (K) [A]
Chole = 4.3e-3;                 % conduct. of vent. holes (m^3/s) [A]
Chy = 1.33;                     % ratio of conductance (y) [A]
Chz = 1.4;                      % ratio of conductance (z) [1.4]
Vp = 0.1;                       % patch voltage (V) [A]
% patch V fluctuation [W2]
dVp = sqrt( (290e-6.*(1e-4./f)).^2 + (80e-6*ones(size(f))).^2 );

% housing design parameters
alphah = 5e-6;                  % CTE of housing (1/K) [A]
Ch = (s/0.046)^2*(0.004/d)*20.99e-12;% TM-housing capac.(LISA scaled) (F)
Cx = (s/0.046)^2*(0.004/d)*2.31e-12; % TM-x elect. capac.(LISA scaled) (F)
LA = 1e-6;                      % TM-hous. diele. loss. ang. (rad) [A]
Mh = 1;                         % mass of housing (kg)
Vd = 0.1;                       % V diff. across oposite faces (V)
V0g = 0.1;                      % ave. elect. volt. wrt gnd. (V) [S]
Xnr = 5e-9;                     % TM sensing error for D-F control [S]
Va = 0;                         % maximum applied voltage
dVa = 2e-6*ones(size(f));       % relative applied voltage noise (d-less)

% IMS parameters
Pl = 1e-4;                      % laser power on test mass (W) [A]
RIN = 1e-3;                     % laser rel. int. noise (1/Hz^1/2) [A]

% spacecraft parameters
M = 500;                        % mass of SC (kg)
A = 0.5*0.3;                    % area of SC (m^2)
dx = 10e-6;                     % SC jitter (m/Hz^1/2)
Mdis = 1;                       % mass moved by CTE (kg) [S]
xdis = 0.3;                     % dist. to mass moved by CTE (m)[SmallSat]

B = 20*5e-6;                    % SC dc magnetic field (T) [10*A]
dB = 20*100e-9;                 % SC mag fld fluct. (T/Hz^1/2) [10*A]
delB = 20*1e-6;                 % SC dc mag field grad (T/m) [10*A]
ddelB = 20*25e-9;               % SC fld.grad.fluct.(T/m/Hz^1/2)[10*A]
dBac = 20*1.0e-8;               % SC max AC field fluct. (T) [A]

% B = 1.5e-6;                     % SC dc magnetic field (T) [R, req]
% dB = 150e-9;                    % SC mag fld fluct. (T/Hz^1/2) [R, req]
% delB = 5e-6;                    % SC dc mag field grad (T/m) [R, req]
% ddelB = 240e-9;                 % SC fld.grad.fluct.(T/m/Hz^1/2)[R,req]
% dBac = 15e-9;                   % SC max AC field fluct. (T) [R, req]

% B = 0.144e-6;                   % SC dc magnetic field (T) [R, pred]
% dB = 21e-9;                     % SC mag fld fluct. (T/Hz^1/2) [R,pred]
% delB = 11.5e-6;                 % SC dc mag field grad (T/m) [R, pred]
% ddelB = 39e-9;                  % SC fld.grad.fluct.(T/m/Hz^1/2)[Rpred]
% dBac = 0.3e-9;                  % SC max AC field fluct. (T) [R, pred]

delsqB = 10*0.2;                % 2nd of SC mag field (T/m^2) [10*A]
Bac = 10*0.5e-6;                % SC max AC field (T) [10*A]
zeta = 1e-2;                    % SC mag shield factor [D]
a0x = 9e-8;                     % SC dc diff. acceleration (m/s^2) [tors pend x3]
dgth = 1e-8;                    % grav. noise coeff. (m/s^2) [A]
alphadis = 2.4e-5;              % CTE relavant to mvd mass (1/K) [Al]
T = 293;                        % temperature (K) [D]
dTd = 0.005*(1e-3./f).^(1/3);    % temp. diff. fluct. (K/Hz^1/2) 
dTsc = 1*(1e-3./f).^(1/3);      % temp. fluct. of SC (K/Hz^1/2) 
P = 10e-6;                       % pressure inside housing (Pa)
Dt = 10e-6;                     % thruster noise (N) [MIT Mips]

% space environment parameters
Bip = 5e-7;                         % interplanetary mag field (T) [D]
dBip = max(4e-8*(1e-3./f), 4e-8);   % IP mag field fluct. (nT/Hz^1/2) [D]
qdot = 260;                         % charging rate (e/sec) [D]
% qdot = 1000;                        % charging rate (e/sec) [R]
% dq = 5.85e-16*(1e-3./f);            % TM charge fluct.  (C/Hz^1/2) [D]
dq = e/pi*sqrt(qdot/2)./f;          % TM charge fluct.  (C/Hz^1/2) [W]
Ecr = 3.2e-11;                      % energy of cosmic ray, proton (J) [S]
ncr = 30;                           % cosmic ray impact rate (Hz) [S]
mcr = 1.7e-27;                      % mass of cos.ray, proton (kg) [S]
mo = 30;                            % 
vorb = 29.78e3;                     % SC orbit velocity (m/sec) [D]


%%%%%%%%%%%%%%%   END PARAMETER INPUTS   %%%%%%%%%%%%%%%

% derived TM paramters
rho = rho*ones(Nf, 1);              % vectorized density
s = s*ones(Nf, 1);                  % vectorized length
d = d*ones(Nf, 1);                  % vectorized gap
chi0 = chi0*ones(Nf, 1);            % vectorized susceptibility
a = s.^2;                           % vectorized test mass x-area (m^2)
vol = s.^3;                         % test mass volume (m^3)
m = vol.*rho;                       % test mass mass (kg)
sigma0 = sigma0*ones(Nf, 1);        % vectorized elec. conductance
% frequency dependent susceptibility [A]
chi = chi0 - i*(s.^2).*mu0.*sigma0.*pi.*f/12;   % cube

% maximum control acceleration
aMaxControl = 1/2/m(1)*Cx/d(1)*Va^2;

% derived IMS parameters
RIN = RIN*ones(size(f));            % vecotrized RIN

% derived spacecraft parameters
dB = dB*ones(size(f));              % vecotrized mag field fluct.
delB = delB*ones(size(f));          % vectorized dc mag field grad 
ddelB = ddelB*ones(size(f));        % vectorized mag fld. grad. fluct.
delsqB = delsqB*ones(size(f));      % vectorized 2nd of SC mag field
Bac = Bac*ones(size(f));            % vectorized SC max AC field

%%%%%   TM-to-SC stiffness   %%%%%

% gravity gradient stiffness
% kgg = 2*G*Mdis/xdis^3;            % [S]
Kgg = 1e-7;                         % conservative value [A]

% image charges
Ks1 = 1./m./Ch./d.^2.*(Cx./Ch)*q^2; % [S]
% ks1 = 1/m/Ch/d*q^2;               % [A]

% charge x Vpatch
Ks2 = 2./m./d.^2.*(Cx./Ch)*q*Vp;    % [S] = [A]

% applied voltages [S]
Ks3 = Cx./m./d.^2.*(1/6*Vd^2 + (Ch/6/Cx)^2*V0g^2);

% patch fields
% gamma = 6;                                % [S]
% ks4 = gamma*Cx*(Cx/Ch)^2/m/d^2*Vp^2;      % [S]
Ks4 = Cx./m./d.^2*Vp^2;                     % [A]

% magnetic stiffness [A]
Kmag = abs(chi)*3*sqrt(3)./rho/mu0.*delB.^2 ...
    + 3*sqrt(3)./m.*delsqB.*(abs(chi).*vol/mu0*B + Mr);

% total stiffness
K = Kgg + Ks1 + Ks2 + Ks3 + Ks4 + Kmag;

% spacecraft jitter due to thruster noise
dxt = Dt/M./(2*pi*f).^2;                    % (m/Hz^1/2)

% TM acceleration noise due to stiffness
% facotrs of 1.5 applied for cross-coupling and readout noise [A]
A0 = 1.5*(0 + 1.5*Xnr)*K;

%%%%%   magnetic disturbances (m/sec^2/Hz^1/2)   %%%%

% SC magnetic field fluctuations [A]
A1 = sqrt(3)*abs(chi)./rho/mu0.*delB.*dB;

% SC magnetic field grad. fluct. [A]
A1a = sqrt(3)./m.*(Mr + chi0.*vol.*B/mu0).*ddelB;

% acmag noise downconverted [A]
A1b = 5*abs(chi)./s./rho/mu0.*Bac.*dBac;

% interplanetary mag field fluct. SC grad [A]
A2 = zeta*sqrt(3)*abs(chi)./rho/mu0.*delB.*dBip;

% interplanetary mag field x charge fluct. [D]
A3 = zeta*dq./m*vorb*Bip;

% interplanetary mag field fluct. lorentz [A]
A4 = zeta*q./m*vorb.*dBip;

% combined interplanetary magnetic forces (linear sum)
daIPmag = A2 + A3 + A4;

% combined SC magnetic forces (linear sum)
daSCmag = A1 + A1a + A1b;

%%%%%   thermal disturbances   %%%%%

% radiometric [A]
A8 = (a/2)*P/2./m/T.*dTd;

% radiation pressure asymmetry [A]
A9 = 8/3*sigma.*(a/2)*T^3.*dTd./m/c;

% asymmetric outgassing [A]
A9a = s.^2.*a*Io*Thetao./m/T^2/Chole/(1 + 2*(Chy + Chz)).*dTd;

% thermal distortion [A]
A9b = a0x*alphah*dTd;

% self grav. dc force fluct. due to housing CTE [A]
A9c = dgth*alphah*dTd;

% self grav. dc force fluct. due to SC CTE [S]
A10 = 2*G*Mdis./xdis.^2*alphadis.*dTsc;

% combined thermal disturbances (linear sum)
daT = A8 + A9 + A9a + A9b + A9c + A10;

%%%%%   cosmic ray impacts   %%%%%
A5 = sqrt(2*mcr*Ecr*ncr)./m;

%%%%%    laser photons [A]   %%%%
dPl = Pl*sqrt(RIN);
A7 = 2*dPl./m/c;

%%%%%   brownian noise   %%%%%

% residual gas impacts computed as:
betaRatio = (log(46/4)*(4/46)^2) ...
    / (log(s(1)/d(1))*(d(1)/s(1))^2);
A6 = 1.3e-15*ones(size(f(:, 1)))*betaRatio*sqrt(P/1e-6);        % [W]

% dielectric losses [A]
A6a = 2*Ch./m./d.*(Vp + q./Ch).*sqrt(4*T*kB*LA/2/pi./f./Ch);

% magnetic/eddy current damping [A]
A6b = (3/2/pi)^(1/3).*sqrt(kB*T*s.^5.*sigma0/5).*delB...
    ./ m.*ones(size(f));

% mag. impurities/viscosity in ferromag. component [A]
A6c = 2*delB/pi./m.*sqrt(kB*T.*s.^3.*dchi/mu0./f);

%%%%%   electric force noise   %%%%%

% Vpatch x Vpatch fluct
% A11 = 2./m./d.*Cx*Vp.*dVp;              % [S]
A11 = 1./m./d.*Cx*Vd.*dVp;              % [R]

% Vpatch fluct. x TM charge
% A12 = 2./m./d.*Cx./Ch*q.*dVp;           % [S]
A12 = 1./m./d.*Cx./Ch*q.*dVp;           % [R]

% TM charge fluct. x Vpatch
% A13 = 2./m./d.*Cx./Ch*Vp.*dq;           % [S]
A13 = 1./m./d.*Cx./Ch*Vd.*dq;           % [R]

% TM charge fluct. x charge
% A14 = 2./m./d.*Cx./Ch.^2*dd./d*q.*dq;   % [S]
A14 = 1./m./d.*Cx./Ch.^2*q.*dq;   % [S]

% actuation cross-coupling disturbance
% Aact = 18.7e-16*ones(size(f));          % [A]
Aact = a0x*dVa;

% combined total disturbance (m/sec^2/Hz^1/2)
da = sqrt( A0.^2 ...                                % stiffness
    + daIPmag.^2 + daSCmag.^2 ...                   % magnetic
    + daT.^2 ...                                    % thermal
    + A5.^2 + A7.^2 ...                             % other
    + A6.^2 + A6a.^2 + A6b.^2 + A6c.^2 ...          % brownian
    + A11.^2 + A12.^2 + A13.^2 + A14.^2 ...         % electrical
    + Aact.^2 );                                    % actuation

% da LISA requirement (m/sec^2/Hz^1/2)
daGRACE3 = 3e-13*sqrt(1 + (0.010./f).^(2/3));


% plot results
figure
loglog(f, da, 'k', 'linewidth', 3)
hold on
% grid on
loglog(f, daGRACE3, 'k--', 'linewidth', 3)
loglog(f, A0, 'k', 'linewidth', 2)
loglog(f, sqrt(daIPmag.^2 + daSCmag.^2), 'm', 'linewidth', 2)
loglog(f, daT, 'b', 'linewidth', 2)
loglog(f, sqrt(A5.^2 + A7.^2), 'g', 'linewidth', 2)
loglog(f, sqrt(A6.^2 + A6a.^2 + A6b.^2 + A6c.^2), 'r--', 'linewidth', 2)
loglog(f, sqrt(A11.^2 + A12.^2 + A13.^2 + A14.^2), 'm--', 'linewidth', 2)
loglog(f, Aact, 'b--', 'linewidth', 2)
xlabel('frequency (Hz)')
ylabel('acceleration noise (m/s^2 Hz^1^/^2)')
legend('total estimate', 'spec performance', 'stiffness', 'magnetic', ...
    'thermal', 'sensing, rays', 'Brownian', 'electrical', 'actuation', ...
    'Location', 'southwest');
% axis([1e-4 1 1e-16 1e-12])


mEH = (3 + 2*0.1 + 1)^2*1*10.28*6;


