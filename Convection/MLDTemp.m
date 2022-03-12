function [T] = MLDTemp(Tconst, Tini, Tfin)
T = (Tfin-Tini)/log((Tconst-Tini)/(Tconst-Tfin));
end