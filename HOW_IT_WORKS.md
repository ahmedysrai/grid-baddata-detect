# How it works

## The problem

A control centre estimates the voltage at every bus from a few hundred
noisy meter readings, and checks the fit to catch meters that are wrong.
The usual checks look at the residual, the part of the readings the
estimate cannot explain. Liu, Ning and Reiter (ACM CCS 2009,
doi:10.1145/1653662.1653666) showed that an attacker who knows the network
model can change a chosen set of meters so that the residual does not move
at all. This repo reproduces that on the IEEE 14-bus case, measures how
well the two standard tests do against ordinary gross errors and against
the stealthy attack, and works out how few protected PMUs close the gap.

## Data flow

```
pandapower case14 power flow (the true state)
  -> 82 SCADA meters: |V|, P, Q at every bus, P, Q at one end of each branch
  -> Gaussian noise, one draw per trial           (gridbad.measurements)
  -> attack or gross error added to chosen meters (gridbad.attack)
  -> WLS estimate x_hat, residual r, J, H         (gridbad.wls)
  -> chi-square test on J, largest normalised residual test (gridbad.chi2,
     gridbad.residuals)
  -> flagged or not, counted over many draws      (experiments/)
```

The state is the 14 voltage magnitudes and 13 angles (the slack angle is
the reference), so n = 27. With m = 82 meters the residual has
m - n = 55 degrees of freedom.

## Methods, and why these ones

Estimation is pandapower's WLS, which minimises

$$J(x) = \sum_i \frac{(z_i - h_i(x))^2}{\sigma_i^2}.$$

gridbad.wls wraps it and pulls out the residual, the Jacobian H at the
estimate and h(x) itself, which the attack code needs. Using pandapower's
solver rather than writing one means the estimator is one many people
already use and check; a test checks that the noise-free estimate gives back the
power flow state (to 1e-8 pu; it actually comes out near 1e-11).

Two bad data tests, both standard:

- Chi-square: flag when J exceeds the 95% point of $\chi^2_{55}$ (73.3).
- Largest normalised residual (LNR): with $S = I - H G^{-1} H^T R^{-1}$
  and $G = H^T R^{-1} H$, the normalised residual is
  $r_i / (\sigma_i \sqrt{S_{ii}})$. Flag when the largest exceeds a
  threshold, and remove that meter and re-estimate to identify it.

LNR is there because J is a single number for all 82 meters and dilutes
an error on a meter the estimate absorbs. The bus 4 P injection has the
lowest $S_{ii}$ of the set, 0.282: only 28% of an error on it stays in its
own residual. At a threshold of 3, LNR fires on about 0.19 of clean draws
(82 residuals), so the comparison in detection_rate.py first calibrates
both thresholds to 5% false alarms on 1000 clean draws. That gave 74.1 for
J and 3.45 for LNR, close to the table value and to the Sidak value of
3.42 for 82 independent residuals.

The attack is the exact AC version, $a = h(\hat x + c) - h(\hat x)$,
alongside the linear $a = Hc$ from the paper. The linear one is only
stealthy to first order on an AC estimator; the exact one leaves the
residual unchanged at any size.

The defence is protected PMUs: meters the attacker cannot write. The
attacker can still use any $c$ with $H_{pmu} c = 0$, so a state is safe
only if no vector in the null space of $H_{pmu}$ moves it
(gridbad.placement). Bobba et al. (Workshop on Secure Control Systems,
2010) show for the DC model that protecting a basic measurement set is
necessary and sufficient; the null space test is the same idea on the AC
Jacobian at one operating point. Placements are found by brute force,
which is cheap on 14 buses (under 10 s for the whole experiment) and gives
every optimal placement rather than one.

A zero injection bus adds two constraints the attacker cannot touch:
P = Q = 0 there whatever the meters say. Bus 7 is the only one in case14.

## What is simulated, and where that flatters the method

The whole project is simulation based. It uses one operating point of
case14, with the case's standard loads.

- Meter noise is independent Gaussian with fixed standard deviations:
  0.004 pu for |V| and 1 MW or 1 Mvar for powers. These are assumptions
  with no source yet. Real SCADA errors include bias, time skew between
  meters and quantisation, none of which is modelled, and all of which
  would raise the false alarm rate at a fixed threshold.
- PMU angle noise is 0.2 deg, taken as about a third of the 0.57 deg phase
  error that uses up the 1% TVE limit in IEEE C37.118.1-2011. The 3 sd
  reading is mine.
- The measurement set is redundant everywhere (every bus has V, P and Q).
  Real systems have far fewer meters and some critical ones, where no
  residual test sees anything; the m = n test case shows that limit, the
  experiments do not explore it.
- The attacker knows H exactly and writes meters with no delay. That is
  the strong attacker from the paper, so the stealth results are a worst
  case for the operator, while the detection results for gross errors are
  a best case.
- Current PMUs are checked only to first order. pandapower's estimator has
  no current angle measurement, so the PMU placement result is a rank
  result at the power flow state, not an AC attack run through the
  estimator.

## Reproducing the numbers

```
pip install -e .[dev]
pytest                                        # 36 tests, about 25 s
python experiments/gross_error.py             # about 5 s
python experiments/detection_rate.py          # about 17 min
python experiments/stealthy_attack.py         # about 4 min
python experiments/pmu_check.py               # about 3 min
python experiments/pmu_placement.py           # about 10 s
python experiments/zero_injection_attack.py   # about 3 min
```

Each script prints the tables quoted in NOTES.md and writes its figure to
figures/. Each has a fixed default seed, which --seed changes.

## Results so far

Gross errors, at matched 5% false alarm rates, 200 noise draws
(detection_rate.py): LNR reaches a 0.5 detection rate about 1.2 sd sooner
than chi-square on a line flow and a |V| meter, and about 2 sd sooner on
the bus 4 injection. At 8 sd on bus 4, LNR flags 0.81 of draws and
chi-square 0.47.

Stealthy attack on the bus 14 angle, 50 draws (stealthy_attack.py): it
touches 10 meters. The exact AC attack is flagged at the clean rate (0.02)
at every size up to 30 deg, where 10 deg alone moves the bus 14 P
injection by about 105 MW against a real load of 14.9 MW. The linear
$a = Hc$ is caught every time from 10 deg. A random error on the same 10
meters with the same weighted size is caught at 1 deg.

One protected PMU angle (pmu_check.py): on bus 14 it catches the attack in
every draw from 2 deg; on neighbouring bus 9 it never does, because the
attack does not change h for that angle.

Placement (pmu_placement.py): angle-only and phasor PMUs pin one and two
states each, so they would need one per bus. PMUs that also report branch
currents pin their bus and its neighbours. Four of them pin all 27 states
(five placements, for example IEEE buses {2, 6, 7, 9}). With the bus 7
zero injection three are enough, and {2, 6, 9} is the only such set.

The bus 8 direction that {2, 6, 9} leaves open without the constraint,
run through the AC estimator (zero_injection_attack.py, 50 draws): with
the bus 7 zero injection both tests catch every draw from 1 deg or
0.02 pu. Smaller shifts get through; 0.25 deg moves the bus 7 injection by
about 3 MW. The PMUs are not in z in that run, and the estimator spreads
the mismatch by moving buses 7 and 9, so I think these rates are a lower
bound.

## Limits and next steps

- The SCADA noise levels need a source before the detection rates mean
  anything outside this setup.
- One operating point. The placement result could change with load, since
  H does; repeating it over a set of load scalings is the obvious check.
- Current PMUs in the estimator need current angle rows in h(x), which
  pandapower does not have. With them the {2, 6, 9} result could be run as
  an AC attack instead of a rank test.
- lnr_identify does not take the zero injection option yet.
- Only case14. The brute force placement search grows quickly with bus
  count; a larger case would need an integer program instead.
