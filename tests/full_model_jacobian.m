function J = full_model_jacobian(~, X, model)
% FULL_MODEL_JACOBIAN Test helper: analytic sparse Jacobian for ode15s.
x = X(1:model.size);
Ef = X(end-1);
J = [model.A0+Ef*model.AE, sparse(model.AE*x), sparse(model.initiation); ...
     model.E_release0+Ef*model.E_release1, sparse(model.E_release1*x), sparse(1, 1); ...
     model.Pol_release0+Ef*model.Pol_release1, sparse(model.Pol_release1*x), sparse(-model.k_in)];
end
