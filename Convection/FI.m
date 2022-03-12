function [fi, n, p] = FI(escalfament, gas, prop_fluid, prop_paret)
if gas == 1
    Tm = prop_fluid;
    Tparet = prop_paret;
    fi = Tm/Tparet;
    if escalfament == 1
        p = 0.52;
        n = 0.47;
    else
        p = 0.38;
        n = 0.36;
    end
else
    mu = prop_fluid;
    mu_paret = prop_paret;
    fi = mu/mu_paret;
    if escalfament == 1
        p = -0.33;
        n = 0.11;
    else
        p = -0.24;
        n = 0.25;
    end
end