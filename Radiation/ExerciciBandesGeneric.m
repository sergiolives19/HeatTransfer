sigma = 5.67e-8;

% Dades
n=4;
numk = 2;
h3 = 80;

% Temperature
Tint = 100 + 273.15;
Text = 25 + 273.15;
T1 = 900 + 273.15;
T2 = 7 + 273.15;
T3 = 820.9;     %suposar valor
T4 = Text;
T = [T1, T2, T3, T4];

% Epsilons
epsilon = [];
epsilon(1,1) = 1;
epsilon(1,2) = 1;
epsilon(2,1) = 0.3;
epsilon(2,2) = 0.85;
epsilon(3,1) = 0.9;
epsilon(3,2) = 0.9;
epsilon(4,1) = 1;
epsilon(4,2) = 1;

A = [];
A(1) = 0.19635;
A(2) = 0.19635;
A(3) = 0.11781;
A(4) = 0.03927;
%%% Factors de Visió
F = [];
F(1,1) = 0;
F(2,2) = 0;
F(1,2) = 0.67208;
F(3,3) = 0.12612;
F(3,4) = 0.05408;
F(4,1) = 0.4099;

F(2,1) = RS(F(1,2), A(1), A(2));
F(1,4) = RS(F(4,1), A(4), A(1));

F(1,3) = 1-F(1,1)-F(1,2)-F(1,4);
F(3,1) = RS(F(1,3), A(1), A(3));
F(3,2) = 1-F(3,1)-F(3,3)-F(3,4);
F(2,3) = RS(F(3,2), A(3), A(2));
F(2,4) = 1-F(2,1)-F(2,3)-F(2,4);
F(4,2) = RS(F(2,4), A(2), A(4));
F(4,3) = RS(F(3,4), A(3), A(4));
F(4,4) = 1-F(4,1)-F(4,2)-F(4,3);

assert(ComprovaF(F)==1)

lambda1 = 1e-20;
lambda2 = 5e-6;
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


ro=zeros(n,1);
for i=1:n
    for k=1:numk
        ro(i,k) = 1 - epsilon(i,k);
    end
end

VectorMatriuCoef = [];
b=zeros(n,numk);
J=zeros(n,numk);
G=zeros(n,numk);
fluxQ = zeros(n,numk);
q=zeros(n,numk);

for k=1:numk
    Mk = [];
    for i=1:n
        fila = zeros(1, n);
        for j=1:n
            if i==j
                fila(j)=1-ro(i,k)*F(i,j);
            else
                fila(j)= -ro(i,k)*F(i,j);
            end
        end
    Mk = vertcat(Mk, fila);
    end
    VectorMatriuCoef(:, :, k) = Mk;
    
    for i=1:n
        b(i,k)= E(i,k);
    end

    J(:,k)=VectorMatriuCoef(:,:,k)\b(:,k);
    
    G(:,k) = F*J(:,k);
    fluxQ(:,k) = G(:,k)-J(:,k);
    q(:,k) = fluxQ(:,k).*A';
end

fluxQconv = h3*(T3-Tint);

error = (fluxQ(3,1)+fluxQ(3,2))-fluxQconv;

E = E(:,1)+E(:,2);
J = J(:,1)+J(:,2);
G = G(:,1)+G(:,2);
fluxQ = fluxQ(:,1)+fluxQ(:,2);
q = q(:,1)+q(:,2);

qtotal = 0;
for i = 1:n
    qtotal = qtotal+q(i);
end
qtotal;
    
error = fluxQ(3) - h3*(T(3)-Tint)