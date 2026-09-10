function dXdt = ode_dynamics_full_multipleE(~, X, model)
% ODE_DYNAMICS_FULL_MULTIPLEE Simultaneous binding, recognition and transport.
% State X = [all R microstates; all RHE microstates; E_free; Pol_free].
% Model is prepared once with build_full_rate_matrices. No solver is nested
% in the RHS; free pools evolve by exact resource bookkeeping.
x = X(1:model.size);
Ef = X(end-1);
Pol_f = X(end);
dx = model.A0*x + Ef*(model.AE*x) + Pol_f*model.initiation;
dEf = model.E_release0*x + Ef*(model.E_release1*x);
dPol = model.Pol_release0*x + Ef*(model.Pol_release1*x) - model.k_in*Pol_f;
dXdt = [dx; dEf; dPol];
end
