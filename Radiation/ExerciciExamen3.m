A = [0.05, 7.01858, 18.84956, 7.06858];

% Factor de visio
F = ones(4,4);
F(1,1) = 0;
F(2,2) = 0;
F(4,4) = 0;
F(1,2) = 0;
F(2,1) = 0;

r1 = sqrt(A(1)/pi);
r4 = sqrt(A(4)/pi);
D = 2*r4;
L = A(3)/(pi*D);

R1 = r1/L;
R4 = r4/L;
S = 1+(1+R4^2)/(R1^2);
F(1,4) = 0.5*(S-(S^2-4*(r4/r1)^2)^0.5);
F(4,1) = RS(F(1,4),A(1),A(4));
F(1,3) = 1-F(1,4);
F(3,1) = RS(F(1,3),A(1),A(3));

r5 = r4;
R4 = r4/L;
R5 = r5/L;
S = 1+(1+R5^2)/(R4^2);
F45 = 0.5*(S-(S^2-4*(r5/r4)^2)^0.5);
F(4,3) = 1-F45;
F(3,4) = RS(F(4,3),A(4),A(3));
F(4,2) = 1-F(4,1)-F(4,3);
F(2,4) = RS(F(4,2),A(4),A(2));
F(2,3) = 1-F(2,4);
F(3,2) = RS(F(2,3),A(2),A(3));
F(3,3) = 1-F(3,1)-F(3,2)-F(3,4);
assert(ComprovaF(F)==1)

T = [1200, 500, 700, 500];
epsilon = [1 1; 0.1 0.8; 0.9 0.2; 1 1];

n = 4;
numk = 2;
sigma = 5.67e-8;

lambda1 = 1e-20;
lambda2 = 70e-6;
lambda3 = inf;
lambda = [lambda1, lambda2, lambda3];

%Funcions de radiació
deltaF = [];

for i=1:n
    for k=1:numk
        deltaF(i,k) = DeltaFuncioRad(lambda(k), lambda(k+1), T(i));
        E(i,k) = epsilon(i,k)*sigma*T(i)^4*deltaF(i,k);
    end   
end

deltaF(2,2);
E(1,1);

epsilon = epsilon(:,1);

assert(length(epsilon)==n); % Comprovem que hi han tantes epsilon com parets
[sizeFi, sizeFj] = size(F);
assert(sizeFi==n & sizeFj==n); % Comprovem que F siga (nxn)

%%- Inicialitzar vectors
E = zeros(n, 1);
fluxQ = zeros(n, 1);
ro=zeros(n,1);
for i=1:n
    ro(i) = 1 - epsilon(i);
end

%Knowns
knownT=[1, 1, 1, 1];
knownQ=[0, 0, 0, 0];
    
for i=1:n
    if knownT(i)==1
        assert(knownQ(i)==0) %Si sabem T comprobar que no sabem QA
        assert(T(i)~=0) %I que la T que sabem no val zero!
        E(i)=PoderEmissiu(epsilon(i),T(i));
    else
        assert(knownQ(i)==1) %Si no sabem T comprobar que sabem QA
    end
end
assert((sum(knownT)+sum(knownQ))==n); %Comprovem que no esta infra/sobredimensionat

MatriuCoef = [];

for i=1:n
    fila = zeros(1, n);
    if knownT(i) == 1
        
           for j=1:n
            if i==j
                fila(j)=ro(i)*F(i,j)-1;
            else
                fila(j)=ro(i)*F(i,j);
            end
        end
    elseif knownT(i)==0        
        for j=1:n
            if j==i
                fila(j)=F(i,j)-1;
            else
                fila(j)=F(i,j);
            end
        end
    end
    MatriuCoef = vertcat(MatriuCoef, fila);
end

b=zeros(n,1);
for i=1:n
    if knownT(i)==1
        b(i)=-E(i);
    else
        b(i)=fluxQ(i);
    end
end

J=MatriuCoef\b;


G = zeros(n,1);
for i=1:n
    suma=0;    
    for j=1:n
        suma=suma+J(j)*F(i,j);
    end
    G(i)=suma;    
end
q = zeros(n,1);
for i=1:n
    fluxQ(i) = G(i)-J(i);
    E(i) = J(i) - ro(i)*G(i);
    T(i) = (E(i)/(sigma*epsilon(i)))^(1/4);
    q(i) = fluxQ(i)*A(i);
end

qtotal = 0;
for i = 1:n
    qtotal = qtotal+q(i);
end
qtotal;

%Potencia radiant que la superficie A3 rep
%directament de la superficie A1
J(1)*F(1,3)*A(1)
%Potencia termica total aportada a les superficies
%del recinte escalfades externament
q(1)+q(3)
%Potencia radiant absorbida per la superficie A3
G(3)*A(3)
