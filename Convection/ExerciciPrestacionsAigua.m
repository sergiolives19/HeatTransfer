%Dades
Nt = 9;
Tr = 50+273.15;
Tini = 25+273.15;
uini = 1;
L = 4;
De = 85e-3;
Dti = 0.018;
Dto = 0.02;
LambdaParet = 400;
Nf = 8;
tf = 0.3e-3;
Lf = 1.4e-3;
hi = 1400;

escalfament = 1;
gas = 0;

%Geometria
Smull = pi*(De^2)/4-Nt*pi*(Dto^2)/4;
Pmull = pi*De+Nt*pi*Dto;
Pt = Nt*pi*Dto;
Dh = 4*Smull/Pmull;
A0 = Nt*pi*Dto*L;

%Aletes
CC = 'A';
Apr = Nt*(pi*Dti-Nf*tf)*L;
Af = Nt*Nf*(2*Lf+tf)*L;
Pf = 2*L;
Sf = tf*L;
[r, m] = Aletes(CC, Lf, Sf, Pf, LambdaParet, hi);

%Suposicions inicials
Tfin = 40+273.15;
T0 = 30+273.15;

%Flux màssic
densitat_ini = DensitatAigua(Tini);
FluxMassic = uini*densitat_ini*Smull;
G = m/Smull;

niter = 5;
for i=1:niter
    Tm = (Tini+Tfin)/2;
    %Propietats fluid
    Cp = CalorEspAigua(Tm);
    densitat = DensitatAigua(Tm);
    mu = MuAigua(Tm);
    mu_paret = MuAigua(T0);
    lambda = LambdaAigua(Tm);
    [fi, n, p] = FI(escalfament, gas, mu, mu_paret);
    
    %Constants adimensionals
    Re = Reynolds(FluxMassic, Dh, mu, Smull);    
    Cf = CoefFric(Re);
    Pr = Prandtl(mu, Cp, lambda);
    Nu = Nusselt(Cf, Re, Pr, Dh, L, fi, n);
    
    ho = CoefConv(Nu, lambda, Dh, Pr, Pt, Pmull);
    Rt = 1/(hi*(Apr+r*Af))+log(Dto/Dti)/(LambdaParet*2*pi*L*Nt)+1/(ho*A0);
    U0 = 1/(Rt*A0);
    Tfin = Tr-(Tr-Tini)*exp(-((U0*A0)/(FluxMassic*Cp)));
    q = FluxMassic*Cp*(Tfin-Tini);
    MLDT = MLDTemp(Tr, Tini, Tfin);
    Tm = Tr-MLDT;
    T0 = Tm+q/(A0*ho);
end
PerduaCarga = (2*Cf*(fi^p)*L*G^2)/(Dh*densitat);