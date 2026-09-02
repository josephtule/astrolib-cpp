# Offline Cooperative OD Demo

Load `offline_cooperative_od_demo.json`. No external data files are needed.
The scenario defines bodies and instruments, not measurement schedules or estimator jobs.
It starts paused; an offline driver must unpause its stepping configuration.

## Geometry and Dynamics

- Earth uses WGS84 shape and zonal gravity through J6. Moon and Sun use point-mass
  gravity. All three translate in the truth simulation.
- Earth starts tilted 23.43928 degrees about x and uses simple spin about its body
  z axis at 7.292115e-5 rad/s. The JSON expresses both angle and angular velocity
  in degrees/degrees per second, as required by the shared input angle unit.
- Initial Earth position and velocity are zero in the simulation inertial frame.
  This is an Earth-centered origin at the initial instant, not a frame constrained
  to follow Earth. Subtract Earth's state for Earth-relative OD states.
- J2000 TT is a reproducible clock only. Moon/Sun states, Earth tilt phase, and the
  satellite constellation are synthetic, not ephemerides for that date.
- `leo_lagging` and `leo_leading`: circular 7178.137 km semimajor axes (800 km
  equatorial altitude), 51.6 degree inclination, common plane. Initial along-track
  separation is 25 km, with the leading satellite ahead by about 0.19955 degrees.
- Six GPS-like references use 26560 km circular orbits and 55 degree inclination.
  Nodes are spaced by 60 degrees, with true anomaly opposite the node angle. This
  deliberately clusters their initial positions around the +x side of Earth near
  both LEOs, rather than representing a global GPS constellation.

## Short Arc

Use **t = 0 through 120 seconds** for the first experiment, with a proposed
5-second measurement cadence. These timings are driver inputs, not parsed JSON
fields. The existing fixed RK4 setup subdivides each 1-second world call into four.

The initial GPS-like geocentric angular separations from the LEO pair are below
about 48 degrees. Over 120 seconds the LEOs move about 7.2 degrees and the reference
satellites about 1 degree. This leaves substantial margin below the roughly
74-degree separation at which a GPS-to-LEO line first ceases to point outward from
the LEO radius. The selected links therefore have a conservative short-arc
Earth-clearance geometry. This is a design bound, not a new runtime visibility
check or a guarantee for arbitrary durations. The clustered geometry is also not
an observability/conditioning guarantee for OD.

Runtime horizon, occultation, field-of-view and illumination policies are deferred
until after the comparison plan. No new visibility implementation or diagnostic
was added for this scenario revision.

## Intended Measurement Links

| Observer | Instrument | Target | Component standard deviations |
|---|---|---|---|
| `ground_ranging` | `lagging_range` | `leo_lagging` | 10 m |
| `ground_ranging` | `lagging_range_rate` | `leo_lagging` | 0.01 m/s |
| `ground_angles` | `lagging_radec` | `leo_lagging` | 5 arcsec per angle |
| `leo_lagging` | `leader_range` | `leo_leading` | 1 m |
| Each `gps_like_*` | `tracking_range` | either LEO target | 10 m |
| Each `gps_like_*` | `tracking_range_rate` | either LEO target | 0.01 m/s |

The lagging satellite has **range only**, not range rate. GPS-like platforms now
provide scalar geometric range and range rate, not full relative-position vectors.
These are still idealized instantaneous measurements, not pseudoranges with clock
biases, carrier phase, light time, or a GPS receiver implementation.

Instrument names document intent; the loader does not automatically bind targets.
The same GPS-like instrument can generate separate events for each LEO target.
Range and range rate are separate instruments because no combined type exists yet.
Covariances are diagonal in canonical km, seconds and radians (squared for variances).
Reference satellite ephemerides are treated as known in the first pass.

## Recommended First Scope

1. Use known reference-satellite observer histories. Estimate the lagging satellite
   from ground and GPS-like data; independently estimate the leader from direct
   GPS-like data. Reserve the inter-satellite range as a prediction/residual check.
2. Use declared perturbed initial guesses for the first batch/EKF exercises.
   A 120-second ground-angle arc is not guaranteed to yield robust Gauss IOD, and
   range/range-rate cannot be passed to the existing angle-based IOD algorithms.
   Test angle-based IOD separately on suitable geometry rather than requiring it
   for the initial chained experiment.
3. Match truth/estimator force models for a correctness baseline, or explicitly
   report the model mismatch. Current OD supports two-body/zonal dynamics, not the
   full moving Moon/Sun dynamics in this world. E4 is the planned implementation
   point for consistent provider/stage-time estimator dynamics.

Step 5 supplies the satellite measurement path, not an automatic scenario-driven
multi-estimator pipeline. Keep this full scenario as the future fixture rather
than presenting the short baseline as completed cooperative OD.

## Full Chained Experiment Later

Propagate truth once and generate reproducible noisy measurements. Estimate the
lagging satellite, then record its estimated state and covariance at measurement
epochs. Supply that estimated history, never truth, when estimating the leader.

If the leader receives **only** range from the lagging satellite, its full orbital
state is not generally recoverable robustly from this near-coorbital short arc.
A formation/orbit prior, additional angular information, or direct GPS-like
measurements is needed for a useful experiment. Keeping direct GPS-like links to
the leader makes the eventual crosslink an additional constraint, not the sole
source of its position information.

Remaining work includes an offline coordinator, estimated-observer history replay,
consistent force/frame handling, and observer-uncertainty treatment. Covariance
inflation alone does not preserve temporal or cross-estimator correlations,
especially when both estimates share the same reference-satellite data.

## Follow-ups

- Add measurement scheduling, truth/estimate separation, and explicit link routing.
- Check conditioning, initial-guess sensitivity, and covariance consistency, not
  just final error or small residuals.
- Add runtime visibility/availability policies after the comparison plan.
- Add drag and SRP after their models exist, including area, mass, drag and
  reflectivity parameters; then eclipses. No unsupported JSON placeholders exist.
