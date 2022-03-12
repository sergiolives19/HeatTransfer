function [ Temperatura ] = TempFilter( SolList )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    for merda= 1:length(SolList)
        if abs(imag(SolList(merda))) <= 0.0000001
            if real(SolList(merda)) > 0
                Temperatura = SolList(merda);
                break
            end
        end
    end
end