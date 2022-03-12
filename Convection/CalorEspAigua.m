function [Cp] = CalorEspAigua(T)
Cp = 2820+11.82*T-0.03502*(T^2)+(3.599e-5)*T^3;
end