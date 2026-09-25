# Propagation V1 Regression Artifacts

`manifest.json` describes the immutable V1 baseline cases and their accepted
non-runtime settings. The CSV files contain only requested state samples:

- `fixed_rk4_legacy.csv`: legacy `step_world_legacy(...)` RK4 path.
- `fixed_rk4_current.csv`: staged `step_world(...)` RK4 path.
- `adaptive_dopri54.csv`: staged adaptive Dormand-Prince 5(4) path.
- `provider_states.csv`: staged fixed-RK4 provider translation/orientation path.

The V1 diagnostic creates a missing artifact but does not overwrite an existing one.
Runtime is printed for manual comparison and is intentionally excluded from the
acceptance contract. Propagation Core V2 must compare against these V1 artifacts;
it must not replace them with V2 output.
