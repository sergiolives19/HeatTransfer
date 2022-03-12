function [ pass ] = ComprovaF( F )
eps=0.0015;
[sizeFi, sizeFj] = size(F);
assert(sizeFi==sizeFj); % Comprovar que files == columnes
n = sizeFi;
resultats = [];
for i=1:n
    suma=0;
    for j=1:n
        suma=suma+(F(i,j));
    end
    resultats(i) = suma;
end

pass = 1;
for i=1:n
    if abs(resultats(i)-1)>= eps
        pass = 0;
    end
end
end



