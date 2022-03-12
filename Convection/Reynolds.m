function [Re] = Reynolds(m, Dh, Mu, S)
G = m/S;
Re = G*Dh/Mu;
if Re>2500
    turbulent = 1;
else
    turbulent = 0;
end