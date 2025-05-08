    %function main_backbone
        global N PAS N_PAS Ef_ss;
        syms Ef real;
        % ------------ MODEL PARAMETERS ------------
        L_a = 100;
        P.k_in    = 2;
        P.kEon    = 0.00025;
        P.kEoff   = 10;
        P.k_e     = 65/L_a;
        P.k_e2    = 30/L_a;
        P.E_total = 70000;
        P.L_total = 100000;
        P.kHon = 0.0004; %not sure
        P.kHoff = 0.1; %not sure
        P.kc = 0.6; %not sure
        P.Ef = P.E_total;
        P.Pol_total = 70000;
        
        P.kPmin   = 1; %not sure
        P.kPmax   = 8; %not sure

        geneLength_bp = 25000;
        PASposition   = 20000;
        N      = floor(geneLength_bp / L_a);  % total nodes
        PAS    = floor(PASposition   / L_a);  % node index of PAS
        N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

        % Position-dependent phosphorylation
        % kP = @(n) P.kPmin + 0.01*n;

        % ------------ STEADY STATES & ODE SETUP ------------
        EBindingNumber = 2; 
        [r_E_BeforePas, r_k_AfterPas] = compute_steady_states(P, EBindingNumber + 1); 

        % Generate data for plotting
        Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp
        RE_vals = sym(zeros(EBindingNumber, PAS));

        for e = 1:EBindingNumber+1
            for i = 1:length(Kp_vals)
                kP_val = Kp_vals(i);
                RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
            end
        end

        P.RE_val_bind_E = simplify(sum(sym(1:(EBindingNumber))' .* RE_vals(2:end, :), 1));
        disp(r_k_AfterPas(1));
        P.kHon_t_Ef = matlabFunction(r_k_AfterPas(1), 'Vars', {Ef});
        P.kHoff_t = r_k_AfterPas(2);
        P.kc_t    = r_k_AfterPas(3); 
        %disp(RE_vals);
        % ------------ SOLVE THE ODE ------------
        % tspan = [0 10000];
        %X0 = zeros(N + N_PAS, 1);
%         X0 = [300*ones(PAS, 1); zeros(2*N_PAS-1, 1)];
% 
% %         try
% %             [t, X] = ode45(@(tt,xx) ode_backbone(tt,xx,P), tspan, X0);
% %         catch ME
% %             fprintf('Simulation stopped: %s\n', ME.message);
% %             % Here you can handle or rethrow ME if you want
% %             return;
% %         end
%         X = fsolve(@(xx) ode_backbone(xx, P), X0);
% 
%         R_sol   = X(1:N);
%         REH_sol = X(N+1:end);
%         disp(Ef_ss);
% 
%         for e = 1:EBindingNumber+1
%             for i = 1:PAS
%                 RE_vals(e, i) = R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss}));
%             end
%         end
% 
%         % Define position axis l_values: node n=PAS => l_values=0
%         l_values =  (1-PAS):(N-PAS);
% 
%         % Plot in one figure, scaling by 100 if you want base pairs offset
%         figure; hold on;
%         for e = 1:EBindingNumber+1
%             plot((1-PAS):0, RE_vals(e,:), 'LineWidth',2);
%         end
%         plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
%         plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
%         xlabel('Gene position (bp relative to PAS)');
%         ylabel('Pol II');
%         legend('show','Location','best');
%         title('Final snapshot of R and REH along the gene');
    %end
