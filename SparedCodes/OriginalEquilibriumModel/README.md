# Original equilibrium-model helpers

These six files preserve the original equilibrium-model implementation:

- `run_termination_simulation.m`: equilibrium binding coupled to polymerase kinetics.
- `ode_dynamics_multipleE.m`: original aggregate R/REH dynamics.
- `calculate_pas_cleavage_profile.m`: original cleavage CDF and CAD calculation.
- `compute_avg_E_bound_numerical.m`: equilibrium E and Ser2P averages.
- `compute_steady_states_numerical.m`: equilibrium state distributions.
- `build_rate_matrix_numerical.m`: original E-binding rate matrix.

Their contents are unchanged by the move. Active analyses use the full finite-rate
helpers in the project root and must not add this archive to their MATLAB path.

The tests in this archive are separate from the active test suite. They temporarily
add this directory and then restore the original path. To run them explicitly from
the project root:

```matlab
runtests('SparedCodes/OriginalEquilibriumModel/tests/test_pas_cleavage_at_pas.m')
```

For manual historical comparisons, explicitly add this directory alongside the
project root. The archived solver still requires Optimization Toolbox (`fsolve`),
and the root's `default_parameters.m` supplies current parameters rather than a
frozen historical parameter set. Use explicit parameter values for comparisons.
