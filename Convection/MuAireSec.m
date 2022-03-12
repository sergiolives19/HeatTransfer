function [mu] = MuAireSec(P, T)
mu = (2.469 + 0.0536*T + P/8280)*10^-6;
end