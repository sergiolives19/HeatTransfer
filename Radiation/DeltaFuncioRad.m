function [ DeltaF ] = DeltaFuncioRad( lambda1, lambda2, T )
DeltaF = FuncioRad(lambda2, T)-FuncioRad(lambda1, T);
end