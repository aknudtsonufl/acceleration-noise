function [B,Bgrad] = onboard_magnetic(dipoles) 
    X0 = 1;
    mu0 = pi*4e-7; %magnetic field constant
    sig0 = 1;
    Nd = 1;
    Nf = 1000;
    f = logspace(-4, 0, Nf)';
    f = f(:, ones(1, Nd));
    X_tm = 2e-5;%X0-j*s*mu0*sig0*pi*f/12;
    M = 1e-9; %check this number

    for a=1:length(dipoles(:,1)) %Question: information given is magnitude of magnetic moment, equation is looking for vector?
        xa = dipoles(a,1:3);
        ma = dipoles(a,4:6);
        rhat = (x-xa)/norm(x-xa);
        Ba = mu0/4/pi*(3*dot(ma,rhat)*rhat-ma)/norm(x-xa)^3; % Equation from page 91 of https://upcommons.upc.edu/bitstream/handle/2117/95156/TMDA1de1.pdf?sequence=1&isAllowed=y
        B = B+Ba;
        Bgrad_a = -3*B/norm(x-xa);
        Bgrad = Bgrad+Bgrad_a;
    end
end