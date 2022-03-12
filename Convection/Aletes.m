function [r, m] = Aletes(CC, Lf, S, P, lambda, h)
m=sqrt(h*P/(lambda*S));
if CC == 'D'
    r = 1/(m*Lf);
elseif CC == 'B'
    r = tanh(m*Lf)/(m*Lf);
elseif CC == 'A'
    Lc = Lf+S/P;
    r = tanh(m*Lc)/(m*Lc);
end