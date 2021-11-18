function [F,T] = Radiometer_ForceandTorque(A,P,T0,dT)
    dF_RdT = A*P/4/T0; %N/K
    F = [dF_RdT.*dT,0,0];
    T = [0,0,0];
end