# am/finitediff: defects to fix

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

## 5. specific_heat lookup reads out of bounds

`specificHeatCalculator` evaluates `specific_heat_temps[i - 1]` at i = 0
(line 1015). A 2-point table (`specific_heat 2 300 480 6000 480`) gives
NaN temperatures. Fix: start the loop at i = 1.

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

## Effect

On ThinWall with matched spot, substrate and fixed faces, this branch
gives a 100 µm × 60 µm pool reaching 50 µm into the substrate (peak about
2,150 K). The energy-conserving conduction source gives 60 µm × 25 µm
(about 1,800 K) at the same power, and overshoots this branch at 2× power.
The example powers (115 W) were presumably tuned with the excess energy
and will need re-tuning once 1 and 2 are fixed.
