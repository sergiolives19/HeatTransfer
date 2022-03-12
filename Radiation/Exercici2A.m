%%-Constants
sigma = 5.67*10^-8;
%%-Dades
n = 2;
epsilon = [0.3, 1]';
T1 = 820;
T4 = 500;   %Anar provant valors
Tfluid = 290;
hc = 50;
%Factors de visió
F = [[0,1];[1,0]];

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
T(2) = T4;

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

for i=1:n
    fluxQ(i) = G(i)-J(i);
    E(i) = J(i) - ro(i)*G(i);
    T(i) = (E(i)/(sigma*epsilon(i)))^(1/4);
end

error = fluxQ(2)-hc*(T4-Tfluid)
fluxQ(2)