% Define the system of equations as a function
function F = mySystem(X)
    R = R0;          
    RE = RE0;     
    REL = REL0;     
    E_f = E_f0;          
    L_f = L_f0;          
    
    F = zeros(length(X), 1);     % Initialize the output vector for all equations
    
    % Define the equation for l = -1000 (special case)
    F(1) = k_in - k_c * R(1) - kE_on * R(1) * E_f + kE_off * RE(1);
    
    % Define the equations for l = -999 to 1000
    for l = 2:N+1
        if l <= 1001  % l < 0, corresponding to indices 1 to 1001
            k_e = 0.9;
        else  % l >= 0, corresponding to indices 1002 to 2001
            k_e = 0.2;
        end
        
        % Index in MATLAB: l = 1 corresponds to l = -1000, l = 2001 corresponds to l = 1000
        F(l) = k_e * R(l-1) - k_e * R(l) - kE_on * R(l) * E_f + kE_off * RE(l);
        F(length(R0) + l) = k_e * RE(l-1) - k_e * RE(l) + kE_on * R(l) * E_f - kE_off * RE(l) + kL_off * REL(l) - kL_on * RE(l) * L_f;
        F(length(RE0) + l) = k_e2 * REL(l-1) - k_e2 * REL(l) - k_c * REL(l) - kL_off * REL(l) + kL_on * RE(l) * L_f;
    end
    
    % Equations for E_f and L_f
    F(length(R0) + length(RE0) + length(REL0) + 1) = E_f - (E_total - sum(RE + REL));
    F(length(R0) + length(RE0) + length(REL0) + 2) = L_f - (L_total - sum(REL));
end
