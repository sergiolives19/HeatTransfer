sigma = 5.67e-8;
D=12e-3;
s=32e-3;
T1=530;
T2=295;
Tinf=480;
hc=200;
A1=s;
A2=s;
A3=pi*D;

F= [[-1,-1,-1];[-1,-1,-1];[-1,-1,-1]];
F(1,1)=0;
F(2,2)=0;
F(1,3)= 1-sqrt(1-(D/s)^2)+(D/s)*atan(sqrt((s^2-D^2)/(D^2)));
F(1,2)= 1-F(1,3);
F(2,1)=RS(F(1,2), A1, A2);
F(2,3)=1-F(2,1);
F(3,1)=RS(F(1,3), A1, A3);
F(3,2)=RS(F(2,3),A2,A3);
F(3,3)=1-F(3,1)-F(3,2);

%Fent un balanç de potencies a la superficie 1
T3 = (T1^4-(hc*(Tinf-T1)+sigma*(T2^4-T1^4)*F(1,2))/(sigma*F(1,3)))^0.25;

T=[T1, T2, T3];

J=[];
for i=1:3    
    J(i)=sigma*1*T(i)^4;
end

G=[];
for i=1:3
    G(i)=J(1)*F(i, 1)+J(2)*F(i, 2)+J(3)*F(i,3);
%    for j=1:3
%        G(i)=G(i)+J(j)*F(j, i);
%    end
end

E1_25 = 1*sigma*T1^4*DeltaFuncioRad(2e-6, 5e-6, T1);