function [Nu] = Nusselt(Cf, Re, Pr, Dh, L, fi, n)
Num = Cf*(Re-1000)*Pr*(1+(Dh/L)^(2/3))*fi^n;
Den = 2+17.96*Cf^0.5*(Pr^(2/3)-1);
Nu = Num/Den;
end