# Changelog

## 1.0.0

First public release of a Solar System N-body code originally written from scratch as a set of
MATLAB scripts, now reorganised as a documented library with examples and tests.

- Vectorised N-body and test-particle accelerations, and a single propagator shared by every scenario.
- Gravitational parameters updated to JPL DE440 values; mean radii for every body.
- Adaptive Dormand-Prince integrator (`ode45`) alongside the original explicit midpoint scheme.
- Close-approach screening generalised to any set of bodies, with miss distance, time of closest
  approach and sub-step resolution of encounters.
- Impact probabilities with binomial uncertainties; calendar dates with leap years.
- Planet positions at the half step interpolated when propagating asteroids.
- Verification and validation suite, including comparison with JPL Horizons, run on GitHub with GNU Octave.
