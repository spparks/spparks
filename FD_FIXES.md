# am/finitediff: defects to fix

**Status (2026-10-16): items 1–7 fixed in the commit after this note** —
energy audit now gives 1.032× the laser energy, flat in time (the +3.2% is
the source counting its centre plane with full weight, see 8); stable with
`k_powder` = `k_solid`; 2-point cp tables work. Laser inputs re-tuned, see
"Recalibration" below.

Found 2026-10-16 while comparing this branch (f6a84fe0) with the conduction
temperature source on SURGE-kokkos (`docs/CONDUCTION_SOURCE.md` there).
Runs and audit decks: `/home/tron/fd_runs` (`ThinWall`, `LargerRectangle`,
`audit`, `cmp`). Line numbers are `src/app_am_finitediff.cpp`.

Built with `cd src && make mpi` in the spparks-dev spack env.

## 1. The thermal update does not conserve energy (most important)

`site_event_finitedifference` writes `T[i] += time_step * chm_temp` in place
(line 980) inside the site loop of `app_update` (line 741), so later sites
read neighbours already advanced this step. The exchange between sites i
and s no longer cancels: it is off by r·ΔT_i, laser heating included
(r = α dt/dx²).

Measured (`fd_runs/audit/in.audit`: constant cp 480, latent 0, k_powder
0.3, losses ~0, 60-site track mid-box, 8 ranks, instrumented copy in
`/home/tron/spparks-fd-audit` printing ΣT every 20 thermal steps): ΣT
rises 2.96× the laser energy over steps 20–40, falling to 1.28×
cumulative by step 280. Expected rise is P dt / (ρ cp dx³) = 2.351e5 per
step.

Fix: compute all new temperatures into a second array, then swap; one
`comm->all()` after. The flux between two sites should be one expression
used with opposite signs (finite-volume form), so it cancels exactly.

## 2. Time step above the explicit stability limit

The examples use `time_step` 7e-7 s with k 30, ρ 5706, cp 480, dx 5e-6:
r = α dt/dx² = 0.31, against the 3-D FTCS limit 1/6 (dt ≤ 3.8e-7 s). It
only stays bounded because of the in-place sweep and the low powder k; with
`k_powder` = `k_solid` it is NaN within 20 steps.

Fix: compute dt_max = dx² ρ cp_min / (6 k_max) in `init_app` and error (or
sub-cycle) when `time_step` exceeds it. Fixing 1 without 2 will blow up.

## 3. Surface losses have the wrong sign

Lines 944–947 and 968–971 add `(h (Tb - T_room) + eta sigma (Tb^4 - T_room^4))
/ density / Cp / dx` to `chm_temp`, i.e. convection and radiation heat a
surface hotter than the room. Should be subtracted. The extrapolated
`Tb = (3 T[i] - T[oppo]) / 2` can overshoot; T[i] is the safe choice.

## 4. Simulation clock stalls while the MC step is short

`iterate_rejection`: when dtMC < time_step, the branch runs the MC sweeps,
then `app_update(time_step)` (line 589), but time only advances
`if (dtMC > 1) time += dtMC` (line 581). Thermal steps (and the laser) go
on while `time` stands still, so output times, `run` lengths and anything
compared against wall-clock time are wrong whenever the pool is hot.

Fix: advance `time += time_step` with each thermal step in that branch.

## 5. specific_heat: table length never set, lookup reads out of bounds

The `specific_heat` parser never sets `specific_heat_length = nHeat`, so
the lookup keeps the default length 6 and reads past a shorter table: a
2-point table (`specific_heat 2 300 480 6000 480`) gives NaN temperatures.
`specificHeatCalculator` also evaluates `specific_heat_temps[i - 1]` at
i = 0 (line 1015). Fix: set the length (and check narg), loop from i = 1.

Related: latent heat enters as an apparent cp over [Ts, Tl] (line 1027). A
step that heats a site by more than Tl − Ts skips it; an enthalpy
formulation (state H, T = T(H)) cannot.

## 6. Solid–powder conductance uses the arithmetic mean

Lines 937 and 960 use (k_solid + k_powder)/2 ≈ 15 W/m·K across a
solid–powder face. Two half-cells in series give the harmonic mean
2 k_s k_p / (k_s + k_p) ≈ 0.6 W/m·K, so powder next to solid conducts
about 25× too well.

## 7. Recoat reset ignores remaining liquid

`app_update` (line 728) resets every site to `boundary_temperature` as soon
as max T < Ts. Liquid sites (activeFlag 2) whose solidification is still
pending then see a jump to room temperature. Also require that no site is
liquid.

## 8. Source normalisation counts the centre plane fully (minor, not fixed)

The half-space factor 2 assumes the Gaussian is cut at its centre plane,
but sites on that plane (z = floor(z_meltspot) when the path z is an
integer) get full weight: +3.2% absorbed power at sigma = 12 cells.

## Effect

On ThinWall with matched spot, substrate and fixed faces, this branch
gives a 100 µm × 60 µm pool reaching 50 µm into the substrate (peak about
2,150 K). The energy-conserving conduction source gives 60 µm × 25 µm
(about 1,800 K) at the same power, and overshoots this branch at 2× power.
The example powers (115 W) were presumably tuned with the excess energy
and will need re-tuning once 1 and 2 are fixed.

## Next to investigate: set `time_step` automatically

`time_step` (app_style arg 3) is the coupled step: laser advance
`vel·time_step`, one thermal update (now sub-cycled to `dt_stable`), one
melt/solidification pass. `init_app` only checks
`max_front_vel·time_step <= dx`. Idea: derive it from that front-speed
limit (with a margin) instead of taking it as input, keeping thermal
sub-cycling to `dt_stable` underneath.

Open question: there were reports of odd behaviour near the front-speed
limit. They may have been thermal-solver artefacts (items 1–4: energy
excess, dt above the FTCS limit, stalled clock), not the limit itself.
Re-test near the limit with the fixed solver before trusting either.

Cost note: with `time_step` 7e-7 s and `dt_stable` 3.42e-7 s each step
takes 3 thermal sub-steps (ThinWall 26 s -> 67 s on 8 ranks). Two would
do if sub-steps were sized to the limit rather than 0.9 of it, and the
per-sub-step `comm->all()` could exchange T only (`all_selective`).

## Recalibration after the fixes (2026-10-16)

With 1–7 fixed, the original inputs (115 W, sigma 60 um) give a 60 x 25 um
ThinWall pool, the same as the conduction source. Re-tuned on ThinWall
(`/home/tron/fd_runs/cal`, 13 runs) to the original pool (median 100 um
wide, 60 um melted depth): **150 W absorbed, x_dev = y_dev = 4.2e-5,
z_dev = 7e-5 m**, applied to all four example decks.

| ThinWall | solid above substrate | pool width / depth |
|---|---|---|
| original (115 W, energy excess) | 58,300 | 100 / 60 um |
| fixed, 150 W, 42/42/70 um | 63,662 | 100 / 60 um |

LargerRectangle (first two layers, fixed, 150 W): 570,300 solid sites
(original 600,600), melted depth 75 / 95 um (original 70 / 85).

Run times now follow the real clock: ThinWall 0.0036 -> 0.005 s (4 layers
complete), SingleLayerScan 0.0361 -> 0.037 s (35.6 ms of scan),
LargerRectangle 0.043 -> 0.0433 s (stops just after layer 2's recoat reset
at 43.25 ms; layer 3 is above the 30-site box). GcodeTranslator unchanged.
