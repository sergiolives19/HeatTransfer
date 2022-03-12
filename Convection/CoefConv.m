function [h] = CoefConv( Nu, lambda, Dh, Pr, Pt, Pmull )
h = (Nu*lambda/Dh)*(1-(0.75/(1+Pr))*(1-(Pt/Pmull)));
end

