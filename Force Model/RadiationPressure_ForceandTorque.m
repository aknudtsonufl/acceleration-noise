function [F,T] = RadiationPressure_ForceandTorque(A,T0,sigma,c,dT)
    dF_RPdT = A*8/3*sigma/c*T0^3; %N/K
    F = [dF_RPdT*dT,0,0];
    T = [0,0,0];
end
