# Legacy Propagation Baseline

This directory stores deterministic state histories generated from
`scenarios/propagation_baseline_demo.json`. The files freeze the accepted legacy
world-stepper behavior before Propagation Core V2 replaces its inner loop.

## Cases

- `legacy_rk4_reference.csv`: fixed RK4 with a 0.1-second nominal step.
- `legacy_dopri54_reference.csv`: adaptive Dormand-Prince 5(4) with the baseline
  tolerances configured in `run_propagation_v1_baseline_diag()`.

Both cases cover 0-1000 seconds and sample Earth, Moon and satellite translation and
attitude at 0, 100, 250, 500, 750 and 1000 seconds. Satellite and Moon orbital
invariants are computed relative to Earth.

## Acceptance

Repeated fixed RK4 runs must be bitwise deterministic through the stored sample
fields. The adaptive case is accepted against fixed RK4 with these absolute limits:

- position: `1e-5 km`
- velocity: `1e-8 km/s`
- quaternion: `1e-9`
- angular velocity: `1e-10 rad/s`
- invariant values: `1e-4` in their stored simulation units

Runtime is informational and is not a cross-machine acceptance threshold. The legacy
stepper exposes derivative evaluations, stage builds and provider queries. Portable
allocation counting is marked unavailable; use an external allocation profiler until
the contiguous V2 workspace provides deterministic allocation guarantees.

The canonical V1 migration artifact set now lives in
`data/regression/propagation_v1/`. These earlier files are retained as historical
references. The normal diagnostic does not overwrite committed artifacts.
