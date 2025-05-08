function equations = generate_equations(n)
    % Define symbolic parameters
    syms R0 Kp Kon Koff Ef
    equations = []; % Initialize equations array

    % Generate symbolic variables for R and RE terms
    R = sym('R', [1, n + 1]);  % R0, R1, ..., Rn
    RE = sym('RE', [n, n]);    % R1E, R2E, ..., RnE, R2E2, ..., RnEn
    
    disp(R);
    % Equation for R0
    eq = Kp * R(2) - R(1) == 0;
    equations = [equations; eq];

    % Generate equations for R and RE terms
    for i = 1:n
        % Equation for Ri
        eq_ri = R(max(i - 1, 1)) + Koff * RE(i, 1); % Start with R(i-1) and RE(i,1)
        if i < n
            eq_ri = eq_ri + Kp * R(i + 1); % Add Kp * R(i+1) if valid
        end
        eq_ri = eq_ri - (Kon * Ef * i + Kp + 1) * R(i);
        eq = eq_ri == 0;
        equations = [equations; eq];

        % Equations for RE terms
        for j = 1:i
            eq_reij = Kon * Ef * R(i);
            if i < n
                eq_reij = eq_reij + Kp * RE(i + 1, j); % Add RE(i+1, j) if valid
            end
            if j > 1
                eq_reij = eq_reij + RE(i - 1, j - 1); % Add RE(i-1, j-1) if valid
            end
            eq_reij = eq_reij - (Kon * Ef + Kp + (j - 1) * Koff) * RE(i, j);
            eq = eq_reij == 0;
            equations = [equations; eq];
        end

        % Final condition for RE terms (when j = i)
        eq_final = Kon * Ef * RE(i, i) - i * Koff * RE(i, i) == 0;
        equations = [equations; eq_final];
    end
end
