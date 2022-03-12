%Dades
a = 0.06;
b = 0.05;
t = 0.003;
L = 0.1;
Pini= 200;
Tini = 40+273.15;
uini = 5;
LambdaParet = 60;
Tv = 100+273.15;
hv = 2550;
hlv = 2257;

escalfament = 1;
gas = 1;

%Suposicions inicials
Tfin = 330;
Tparet = 340;
Pfin = (Pini/Tfin)*Tini;
Tm = (Tini+Tfin)/2;
Pm = (Pini+Pfin)/2;

%Propietats fluid
densitat_ini = DensitatAireSec(Pini, Tini);
densitat = DensitatAireSec(Pm, Tm);
Cp = 1004;
lambda = LambdaAireSec(Tm);
mu = MuAireSec(Pm, Tm);

[fi, n, p] = FI(escalfament, gas, Tm, Tparet);

%Geometria
Lf = (b-3*t)/2;
Smull = 4*Lf^2;
Pmull = 16*Lf;
Pt = Pmull;
Dh = 4*Smull/Pmull
A1 = 4*b*L;
Am = (b-t)*4*L;
A2 = 16*Lf*L;

%Flux màssic
FluxMassic = uini*densitat_ini*Smull;

%Constants adimensionals
Re = Reynolds(FluxMassic, Dh, mu, Smull);    
Cf = CoefFric(Re);
Pr = Prandtl(mu, Cp, lambda);
Nu = Nusselt(Cf, Re, Pr, Dh, L, fi, n);

hi = CoefConv(Nu, lambda, Dh, Pr, Pt, Pmull)

%Aletes
CC = 'B';
Apr = 8*Lf*L;
Af = 8*Lf*L;
Pf = 2*L;
Sf = t*L;
[rendiment, m] = Aletes(CC, Lf, Sf, Pf, LambdaParet, hi);

Rt = 1/(hv*A1)+t/(LambdaParet*Am)+1/(hi*(Apr+rendiment*Af));

U1 = 1/(Rt*A1)
U2 = U1*A1/A2;

Tfin = Tv-(Tv-Tini)*exp(-(U2*A2/(FluxMassic*Cp)));

T = MLDTemp(Tv, Tini, Tfin);

qconv1 = FluxMassic*Cp*(Tfin-Tini)
qconv2 = U2*A2*T

Tparet = Tv-qconv1*(1/(hv*A1)+t/(LambdaParet*Am))

niteracions = 3;
for i=1:niteracions
    Pfin = (Pini/Tfin)*Tini;
    Tm = (Tini+Tfin)/2;
    Pm = (Pini+Pfin)/2;
    fi = Tm/Tparet;
    densitat = DensitatAireSec(Pm, Tm);
    Cp = 1004;
    lambda = LambdaAireSec(Tm);
    mu = MuAireSec(Pm, Tm);
    Re = Reynolds(FluxMassic, Dh, mu, Smull);    
    Cf = CoefFric(Re);
    Pr = Prandtl(mu, Cp, lambda);
    Nu = Nusselt(Cf, Re, Pr, Dh, L, fi, n);
    hi = CoefConv(Nu, lambda, Dh, Pr, Pt, Pmull);
    [rendiment, m] = Aletes(CC, Lf, Sf, Pf, LambdaParet, hi);
    Rt = 1/(hv*A1)+t/(LambdaParet*Am)+1/(hi*(Apr+rendiment*Af));
    U1 = 1/(Rt*A1);
    U2 = U1*A1/A2;
    Tfin = Tv-(Tv-Tini)*exp(-(U2*A2/(FluxMassic*Cp)));
    Pfin = (Pini/Tfin)*Tini;
    Tm = (Tini+Tfin)/2;
    Pm = (Pini+Pfin)/2;
    qconv1 = FluxMassic*Cp*(Tfin-Tini);
    Tparet = Tv-qconv1*(1/(hv*A1)+t/(LambdaParet*Am));
end

Tfin
qconv1
Tparet
mcond = (qconv1/hlv)*(3600/1000)