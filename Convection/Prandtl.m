function [Pr] = Prandtl(mu, Cp, lambda)
Pr = mu*Cp/lambda;
end