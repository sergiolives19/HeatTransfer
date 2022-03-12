% Dades
sigma = 5.67*10^-8;
gruix = 0.01;
lambda = 0.2;
epsilon = [0.3, 0.3, 0.3, 1]';

Tfluid = 290;
hc = 50;
T1 = 820;

% Suposem T4 inicial
T4prov = 750;

for iteracions = 1:5
    syms T3
    T2 =(sigma*epsilon(3)*(T3^4-T4prov^4)+hc*(T3-Tfluid))/(lambda/gruix)+T3;
    equacioT3 = sigma*(T1^4-T2^4)/(1/epsilon(1)+1/epsilon(2)-1)+hc*(Tfluid-T2)+lambda*(T3 - T2)/gruix == 0;
    solucioT3 = vpasolve(equacioT3, T3);
    T3prov = TempFilter(solucioT3);
    T2 = (sigma*epsilon(3)*(T3prov^4-T4prov^4)+hc*(T3prov-Tfluid))/(lambda/gruix)+T3prov;
    syms T4
    equacioT4 = sigma*epsilon(3)*(T3prov^4-T4^4)+hc*(Tfluid-T4) == 0;
    solucioT4 = vpasolve(equacioT4, T4);
    T4prov = TempFilter(solucioT4);
end

T2
T3prov
T4prov

qrad2 = sigma*((T1^4-T2^4)/(1/epsilon(1)+1/epsilon(2)-1))
qconv2 = -hc*(Tfluid-T2)
qcond23 = lambda*(T2-T3prov)/gruix
qconv3 = -hc*(Tfluid-T3prov)
qrad4 = sigma*epsilon(3)*(T3prov^4-T4prov^4)
qconv4 = -hc*(Tfluid-T4prov)