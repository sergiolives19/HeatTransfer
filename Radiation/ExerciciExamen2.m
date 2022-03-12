%Dades
T2 = 600;
Tinf = 250;
L = 0.3;
D = 0.2;
s = 0.12;
lambda = 0.22;
gruix = 0.01;

A1 = pi*D^2/4;
A2 = pi*D*L+pi*(D/2)*s;
A3 = A1;

F = [];
F(1,1) = 0;
F(1,2) = 1;
F(2,1) = RS(F(1,2),A1,A2);
F(2,2) = 1-F(2,1);

%Suposició
T1 = 200;

lambda1 = 1e-20;
lambda2 = 200e-6;
lambda3 = inf;

%Funcions de radiació
deltaF = [];
deltaF(1,1) = DeltaFuncioRad(lambda1,lambda2,T1);
deltaF(1,2) = DeltaFuncioRad(lambda2,lambda3,T1);
deltaF(2,1) = DeltaFuncioRad(lambda1,lambda2,T2);
deltaF(2,2) = DeltaFuncioRad(lambda2,lambda3,T2);

epsilon = [0.1,0.78,0.1];

syms T3
T1 = (A3*sigma*epsilon(3)*(T3^4-Tinf^4))*gruix/(lambda*A1)+T3;
equacioT3 = sigma*(T2^4-T1^4)/(1/(epsilon(1)*A1)+1/(epsilon(2)*A2)-1/A2)-(T1 - T3)/(gruix/(lambda*A1))==0;
solucioT3 = vpasolve(equacioT3, T3);
T3prov = TempFilter(solucioT3);
T1 = (A3*sigma*epsilon(3)*(T3prov^4-Tinf^4))*gruix/(lambda*A1)+T3prov;

T1
T3 = T3prov

%%-Constants
sigma = 5.67*10^-8;

%%-Datos
    n=2;
    epsilon = [0.1, 0.78]';

    %Arees
    A = [A1, A2]';

assert(length(epsilon)==n); % Comprovem que hi han tantes epsilon com parets
[sizeFi, sizeFj] = size(F);
assert(sizeFi==n & sizeFj==n); % Comprovem que F siga (nxn)

%%- Inicialitzar vectors
T = zeros(n, 1);
E = zeros(n, 1);
fluxQ = zeros(n, 1);
ro=zeros(n,1);
for i=1:n
    ro(i) = 1 - epsilon(i);
end

    %Knowns
    knownT=[1, 1];
    knownQ=[0, 0];

    %Temperatures
    T(1) = T1;
    T(2) = T2;
    
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
% for i=1:4
%     if T(i) ~= 'n'
%         J(i) = ro(i)*(J(1)*F(i,1)+J(2)*F(i,2)+J(3)*F(i,3)+J(4)*F(i,4))+epsilon(i)*sigma*T(i)^4;
%     else
%         fluxQ(i)=J(1)*F(i,1)+J(2)*F(i,2)+J(3)*F(i,3)+J(4)*F(i,4)-J(i);
%     end       
% end

for i=1:n
    fila = zeros(1, n);
    if knownT(i) == 1
        
        %J(i) = ro(i)*(J(1)*F(i,1)+J(2)*F(i,2)+J(3)*F(i,3)+J(4)*F(i,4))+E(i);
        % J(1)*(ro(1)*F(1,1)-1)+J(2)*ro(1)*F(1,2)+J(3)*ro(1)*F(1,3)+J(4)*ro(1)*F(1,4) = -E(1);
        %%ro(1)*F(1,1)-1, ro(1)*F(1,2), ro(1)*F(1,3), ro(1)*F(1,4)]
        % J(1)*ro(4)*F(4,1)+J(2)*F(4,2)+J(3)*ro(4)*F(4,3)+J(4)*(ro(4)*F(4,4)-1)= -E(4);
        %%[ro(4)*F(4,1), ro(4)*F(4,2), ro(4)*F(4,3), ro(4)*F(4,4)-1]]
        for j=1:n
            if i==j
                fila(j)=ro(i)*F(i,j)-1;
            else
                fila(j)=ro(i)*F(i,j);
            end
        end
    elseif knownT(i)==0        
        % F(2, 1)*J(1) + (F(2, 2)-1)*J(2) + F(2, 3)*J(3) + F(2, 4)*J(4) = fluxQ(2);
        % F(3, 1)*J(1) + F(3, 2)*J(2) + (F(3, 3)-1)*J(3) + F(3, 4)*J(4) = fluxQ(3);
        %%[F(2,1), F(2,2)-1, F(2,3), F(2,4)];
        %%[F(3,1), F(3,2), F(3,3)-1, F(3,4)];
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
    % G(3) = J(1)*F(3,1)+J(2)*F(3,2)+J(3)*F(3,3)+J(4)*F(3,4);
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