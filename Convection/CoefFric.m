function [Cf] = CoefFric(Re)
%if (1E4<Re<1E6)
%    EsCompleixFilonenko = 1
%else
%    EsCompleixFilonenko = 0    
Cf = 1/(1.58*log(Re)-3.28)^2;
end