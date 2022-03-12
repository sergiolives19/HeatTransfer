function [mu] = MuAigua(T)
mu = 10^(-13.73+1830/T+0.0197*T-(1.47e-5)*T^2);
end