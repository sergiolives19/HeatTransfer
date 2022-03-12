%Dades
a = 0.04;
b = 0.03;
t = 0.002;
Lf = 0.006;
LambdaParet = 280;
Pini= 250;
FluxMassic = 0.03;
Tini = 25+273.15;
q = 1200;
Tv = 100+273.15;
hv = 2500;

escalfament = 1;
gas = 1;

%Geometria
Smull = a^2-b^2;
Pmull = 4*a+4*b;
Pt = 4*b;
Dh = 4*Smull/Pmull;    
A1_prima = 4*b;
Am_prima = 4*(b-t);

%Aletes
CC = 'A';
Apr_prima = 4*(b-2*t)-4*t;
Af_prima = (8*Lf+4*t);
Pf_prima = 2;
Sf_prima = t;
[r, m] = Aletes(CC, Lf, Sf_prima, Pf_prima, LambdaParet, hv);

G = FluxMassic/Smull;

%Propietats fluid
Cp = 1004;
Tfin = q/(FluxMassic*Cp)+Tini;
MLDT = MLDTemp(Tv, Tini, Tfin);
Tm = Tv-MLDT;
Pm = (Pini/Tm)*Tini;
densitat = DensitatAireSec(Pm, Tm);
lambda = LambdaAireSec(Tm);
mu = MuAireSec(Pm, Tm);

%Suposicions
L = 10;
T1 = 350;

niter = 10;

for i=1:niter
    %Constants adimensionals
    [fi, n, p] = FI(escalfament, gas, Tm, T1);
    Re = Reynolds(FluxMassic, Dh, mu, Smull);    
    Cf = CoefFric(Re);
    Pr = Prandtl(mu, Cp, lambda);
    Nu = Nusselt(Cf, Re, Pr, Dh, L, fi, n);

    ho = CoefConv(Nu, lambda, Dh, Pr, Pt, Pmull);

    Rt_prima = 1/(ho*A1_prima)+t/(LambdaParet*Am_prima)+1/(hv*(Apr_prima+r*Af_prima));
    L = q*Rt_prima/MLDT;    
    A1 = 4*b*L;
    T1 = Tm+q/(ho*A1);
end
U1 = 1/(Rt_prima*A1_prima);
PerduaCarga = (2*Cf*(fi^p)*L*G^2)/(Dh*densitat);