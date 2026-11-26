# Notes

## 2026-11-18

First end-to-end run: pandapower's WLS on case14, 82 measurements (V, P, Q
at every bus, P, Q at the sending end of all 20 branches), n = 27, so 55
degrees of freedom and a 95% threshold of 73.3.

Noise-free measurements did not give back the power flow state at first:
0.005 pu and 0.23 degrees off with all 82, and 0.012 pu, 0.39 degrees with
only V and the bus injections. Dropping one bus at a time showed it was Q
at bus 9, where case14 has a 19 Mvar capacitor (21.2 Mvar at the 1.056 pu
it sits at). `res_bus.q_mvar`
includes the shunt, but the estimator puts the shunt in Y, so the injection
"measurement" has to leave it out. With that fixed the state matches to
1e-11 pu.

Gross error of k sd on one measurement, same noise (seed 1, clean J = 44.0;
experiments/gross_error.py):

```
P injection, bus 4    5 sd: J =    50.9    10 sd: J =    72.0    20 sd: J =   156.3*
Q flow, line 2-3      5 sd: J =    60.0    10 sd: J =   118.6*   20 sd: J =   363.6*
V magnitude, bus 14   5 sd: J =    64.5    10 sd: J =   126.2*   20 sd: J =   373.3*
```

The bus 4 injection hides much better: a 10 sd error still passes. Bus 4
has five branches (2-4, 3-4, 4-5 and the transformers to 7 and 9), each
with its flow measured at one end. Probably the estimate absorbs
most of an injection error into the bus 4 state instead of leaving it in
the residual.
Worth computing the residual sensitivity matrix S = I - H G^-1 H^T R^-1
next, so this can be read off as S_ii instead of guessed.

Seed 2 starts from J = 60.7 and flags more of the 5 sd cases, so one seed
says little about detection rate at small k. Needs a Monte Carlo sweep.

## 2026-11-19

Added S_ii and normalised residuals (src/gridbad/residuals.py), with H
taken at the estimate. On the full 82-measurement set the S_ii sum to 55
as they should, and range from 0.282 to 0.958.

The bus 4 P injection is the 0.282, the lowest of all 82, so yesterday's
guess holds: only 28% of an error on it stays in its own residual.
Same seed 1 noise, r_N of the measurement that carries the error
(experiments/gross_error.py):

```
                      S_ii   r_N at k = 5, 10, 20
P injection, bus 4   0.282     2.6     5.3*   10.6*
Q flow, line 2-3     0.854     4.0*    8.7*   17.9*
V magnitude, bus 14  0.822     4.5*    9.1*   18.1*
* largest r_N in the set and above 3
```

For a single error with no noise r_N = k sqrt(S_ii): 5.3 for bus 4 at
k = 10, which is what came out. The 10 sd bus 4 error that J missed (72.0
against 73.3) is caught by LNR, and in every case above 3 the flagged
measurement is the one with the error. With seed 2 the 5 sd bus 4 case
also flags (3.6), so small k still needs the Monte Carlo sweep.

lnr_identify drops the same Q flow as pandapower's perform_rn_max_test on
a 20 sd error, then stops with every r_N below 3.

A measurement set with m = n (V everywhere, P at the 13 non-slack buses)
has every S_ii below 2e-14: all critical, no residual test sees anything
(tests/test_residuals.py).

## 2026-11-20

Monte Carlo sweep over k, 200 noise draws (seed 0), same three targets,
common noise across tests and targets (experiments/detection_rate.py,
figures/detection_rate.png). Fraction of trials flagged:

```
                       k    0     4     6     8    10
P injection, bus 4  chi2  0.04  0.09  0.20  0.48  0.75
                    LNR   0.19  0.39  0.72  0.90  0.98
Q flow, line 2-3    chi2  0.04  0.30  0.75  0.97  1.00
                    LNR   0.19  0.77  0.99  1.00  1.00
V magnitude, bus 14 chi2  0.04  0.34  0.81  0.98  1.00
                    LNR   0.19  0.84  1.00  1.00  1.00
```

LNR looks better everywhere, but the k = 0 column says the comparison is
not fair yet. With no bad data at all, max r_N > 3 in 37 of 200 trials.
That is what 82 roughly independent N(0, 1) residuals give:
1 - (1 - 0.0027)^82 = 0.20. The chi-square test is at 0.035 here (0.05
expected, 7 of 200). So LNR at 3 runs at about four times the false alarm
rate, and part of its lead is just that.

On the bus 4 injection the gap is still large where the false alarms
cannot explain it: at k = 8, 0.83 of trials have LNR flag the right
measurement, against 0.48 for chi-square flagging anything. On the other
two the curves are a couple of sd apart.

To compare at the same false alarm rate, the LNR threshold for 5% overall
is about 3.42 if the 82 residuals were independent (per-measurement
two-sided alpha 0.05/82). They are correlated through S, so the real
value is a little lower; it should come from the k = 0 runs instead.

## 2026-11-21

Calibrated both thresholds on 1000 clean noise draws (separate stream from
the test trials) and reran the sweep at matched false alarm rates
(experiments/detection_rate.py, figures/detection_rate.png redrawn):

```
95% point of J(x)     74.1   (chi-square table, 55 dof: 73.3)
95% point of max r_N  3.45   (Sidak, 82 independent residuals: 3.42)
```

I expected the empirical LNR value to come out below 3.42, since the
correlation through S should make the maximum a bit smaller than for
independent residuals. It came out 0.03 above. With 1000 draws the 95%
point is not pinned down that finely, so all this says is that the
correlation effect is small here and the Sidak value is a fair default.
On the 200 test trials the false alarm rates were 0.030 (chi-square) and
0.020 (LNR), both a little under 0.05 but within about two binomial
standard deviations of it.

At matched rates, fraction of trials flagged:

```
                       k    2     4     6     8    10
P injection, bus 4  chi2  0.06  0.08  0.17  0.47  0.72
                    LNR   0.04  0.17  0.47  0.81  0.95
Q flow, line 2-3    chi2  0.06  0.29  0.73  0.96  1.00
                    LNR   0.07  0.56  0.97  1.00  1.00
V magnitude, bus 14 chi2  0.04  0.30  0.79  0.98  1.00
                    LNR   0.06  0.64  0.98  1.00  1.00
```

So the lead LNR had on 2026-11-20 holds after the false alarm correction;
the old numbers at k = 4-6 were inflated, the bus 4 gap is not. Reading the
0.5 crossing off by linear interpolation, LNR gets there about 1.2 sd
earlier than chi-square on the flow and voltage, and about 2 sd earlier
on the bus 4 injection. When LNR fires it nearly
always names the right measurement (0.79 of 0.81 at k = 8 on bus 4).

## 2026-11-22

Stealthy attack on the bus 14 angle (`experiments/stealthy_attack.py`,
50 draws, seed 0). The attacker adds c to the estimated angle of bus 14
and changes every meter that sees it. That is 10 meters: P and Q
injection at buses 9, 13 and 14, and P and Q at the from end of lines
9-14 and 13-14. Nothing else in z depends on that one angle.

The size of what it hides is large. At 10 degrees the exact attack moves
the bus 14 P injection by 105 MW and the flow on 9-14 by 56 MW, against a
real bus 14 load of 14.9 MW.

First draw, thresholds J 73.3 and max r_N 3.42, clean J 58.0 and max r_N
2.41:

```
shift  attack       J   max r_N  angle moved (deg)
    5  linear    74.7      4.67       4.98
    5  ac        58.0      2.41       4.99
    5  random  2843.7     42.16       2.20
   20  linear  2707.2     47.67      19.07
   20  ac        57.9      2.42      19.97
   20  random 50634.2    133.32       1.31
```

Fraction of draws flagged (chi-square / LNR):

```
shift      linear           ac       random
    0  0.02 / 0.02  0.02 / 0.02  0.02 / 0.02
    1  0.02 / 0.02  0.02 / 0.02  1.00 / 1.00
    2  0.02 / 0.02  0.02 / 0.02  1.00 / 1.00
    5  0.26 / 0.26  0.02 / 0.02  1.00 / 1.00
   10  1.00 / 1.00  0.02 / 0.02  1.00 / 1.00
   20  1.00 / 1.00  0.02 / 0.02  1.00 / 1.00
   30  1.00 / 1.00  0.02 / 0.02  1.00 / 1.00
```

The exact AC attack is invisible to both tests at every size tried, up to
30 degrees: J and the largest normalised residual stay at their clean
values to the second decimal and the false alarm rate stays at the clean
0.02. The DC-style a = H c only holds up to a few degrees on the AC model;
the second-order terms it leaves out start to show at 5 degrees and it is
caught every time at 10. The random attack uses the same 10 meters and
the same weighted size as the AC one, and both tests catch it at 1
degree. So stealth comes from matching the column space of H, not from
being small.

On the first try, with 100 draws and a 15 degree step, one attacked
estimate did not converge and the script stopped. I didn't record which
attack it was. A failure now counts as flagged, since the operator would
see it; none happened in the run above.

Next: a detector that does not rely on the residual. A PMU angle at bus
14 or 9 is the obvious one to try, since none of the 10 attacked meters
measures an angle directly.

## 2026-11-23

Added one PMU voltage angle to the 82 SCADA meters and reran the exact AC
attack on the bus 14 angle, with the PMU out of the attacker's reach (its
entry of a set to zero). PMU noise is 0.2 deg: C37.118.1-2011 says a 0.57
deg phase error alone uses up the 1% TVE limit, and I took that as about
3 sd. (experiments/pmu_check.py, 50 draws, about 3 minutes)

First draw, 2 deg shift, PMU at bus 14:

```
 shift  PMU owned  meters      J   max r_N
    2  no            10   144.7     9.31
    2  yes           11    58.6     2.37
```

Fraction of draws flagged (chi-square / LNR):

```
shift          no PMU    PMU at bus 9   PMU at bus 14
    0     0.04 / 0.02     0.04 / 0.04     0.04 / 0.04
 0.25     0.04 / 0.02     0.04 / 0.04     0.06 / 0.02
  0.5     0.04 / 0.02     0.04 / 0.04     0.12 / 0.10
    1     0.04 / 0.02     0.04 / 0.04     0.50 / 0.78
    2     0.04 / 0.02     0.04 / 0.04     1.00 / 1.00
    5     0.04 / 0.04     0.04 / 0.04     1.00 / 1.00
```

A PMU on the attacked bus works: caught every time from 2 deg, and LNR
does better than chi-square at 1 deg (0.78 vs 0.50), as expected when all
the inconsistency sits in one row. If the attacker also owns the PMU it
is one more meter to write (11 instead of 10) and J goes back to the clean
value.

A PMU at bus 9 does nothing, even though bus 9 is next to bus 14. The
attack only moves the bus 14 angle, so h for the bus 9 angle does not
change and the attacker never has to touch that PMU. So a PMU protects
the states it measures, not the area around it. (In the chi-square panel
of figures/pmu_check.png the bus 9 curve sits on top of the no-PMU one.)

The gap below 1 deg is not small in MW terms. On the noise-free estimate
the AC attack changes the bus 14 P injection by 5.1 MW at 0.5 deg and
10.1 MW at 1 deg, against a real load of 14.9 MW, and 0.5 deg is caught
in only about 1 draw in 10. A 1 deg shift is 5 PMU sd, but the
estimator spreads the conflict over the 10 SCADA meters too, so the PMU
residual is smaller than the shift.

Next: which buses need a PMU so that no attack on the SCADA meters alone
can move any state. One PMU per bus is the trivial answer; the question
is how few will do on case14.

## 2026-11-24

How few PMUs stop every SCADA-only attack on case14? The attacker owns
all 82 SCADA meters and none of the PMUs, so it can use any c with
H_pmu c = 0, where H_pmu is the PMU rows of the Jacobian. A state is
safe only if every such c leaves it alone (gridbad/placement.py; Bobba et
al., SCS 2010, show for the DC model that protecting a basic measurement
set is necessary and sufficient). Brute force over every placement, three
PMU models (experiments/pmu_placement.py, about 10 s):

```
PMUs   angle   phasor   currents   best with currents
   1       1        2         12   {4}
   2       2        4         20   {4, 6} {4, 13} {6, 9}
   3       3        6         25   {2, 6, 9}
   4       4        8         27   {2, 6, 7, 9} {2, 6, 8, 9} {2, 7, 10, 13} +2
```

(states pinned out of 27, IEEE bus numbers)

An angle-only PMU, like the one on 2026-11-23, pins one state, and a
voltage phasor two. So with those, protecting everything means a PMU on
every bus. The branch currents are what make placement worth doing: with
the line parameters known, the current phasor and the local voltage give
the voltage at the far end, so each PMU pins its bus and all its
neighbours. Four PMUs pin all 27 states, five placements do it, and no
three do.

The tests check that the pinned set from the null space is exactly the
buses covered by the PMUs and their neighbours, for every placement of
2 and 3, so the rank result and the graph picture agree here.

The best three, {2, 6, 9}, leave only bus 8 open (its angle and |V|).
Bus 8 hangs off bus 7 alone. A null direction scaled to 1 deg there has
|H_pmu c| = 5e-17 and needs 7 SCADA meters written.

This is first order: H_pmu is taken at the power flow state, and
pandapower's estimator has no current angle measurement (the "ia" rows
are commented out in its h(x)), so I could not run the AC attack against
current PMUs the way pmu_check.py does with angles. Bus 7 has no load
or generator; a zero injection constraint there is something the
attacker cannot write either, and might bring the count below 4. Not
tried.

## 2026-11-25

Tried the zero injection idea from yesterday. Bus 7 has nothing on it but
branches (zero_injection_buses finds only bus 7 in case14), so P_7 = Q_7 = 0
is known whatever the bus 7 meters say. An attacker who writes SCADA but
not the PMUs now also has to keep those two at zero, which adds two rows
dS_7/dx to H_pmu before the null space is taken. The rows are checked
against pandapower's own Jacobian rows for the bus 7 P and Q meters.

```
PMUs   currents   currents + bus 7 ZI
   1         12                    14
   2         20                    22
   3         25                    27   {2, 6, 9} only
   4         27
```

So three PMUs with branch currents are enough, and {2, 6, 9} is the only
placement of three that does it. The constraint adds exactly two states
each time, which fits: two equations, and with buses 4, 7 and 9 pinned the
only unknowns left in S_7 are the angle and |V| of bus 8, the bus that
hangs off bus 7 alone. The 1 deg direction that was left open yesterday
moves the bus 7 injection by -11.46 MW and -2.00 Mvar, which the attacker
cannot hide at a bus where it has to be zero.

Same limit as before: first order, at the power flow state, no AC attack
run through the estimator. pandapower can add zero injection as a virtual
measurement (zero_injection= in estimate(), which replaces the bus P and Q
with 0 and a very small std); gridbad.wls passes None so that every row of
z is a real meter. Running the bus 8 attack with that on is the obvious
check.

## 2026-11-26

Ran the bus 8 attack through the AC estimator, with and without the bus 7
zero injection. gridbad.wls.estimate now takes zero_injection=[6]; in
pandapower that does not add rows, it swaps the bus 7 P and Q meters for
virtual zeros at 0.001 pu. So m stays 82 and the degrees of freedom stay
55, which I had wrong yesterday (I expected two extra rows). The attacker
cannot write the virtual rows, so apply_attack skips them.

The attack is the exact AC one on bus 8 only, since with current PMUs on
{2, 6, 9} that is the only direction left (first order). The PMUs are not
in z, because pandapower has no current angle rows. 50 noise draws, seed 0,
thresholds J 73.3 and max r_N 3.42 (experiments/zero_injection_attack.py,
3 minutes):

```
shift          plain (chi2 / LNR)   bus 7 ZI (chi2 / LNR)
none           0.08 / 0.02          0.06 / 0.02
va 0.25 deg    0.08 / 0.02          0.10 / 0.12
va 0.5 deg     0.08 / 0.02          0.38 / 0.66
va 1 deg       0.08 / 0.02          1.00 / 1.00
vm 0.005 pu    0.08 / 0.02          0.10 / 0.22
vm 0.01 pu     0.08 / 0.02          0.54 / 0.86
vm 0.02 pu     0.08 / 0.02          1.00 / 1.00
```

Without the constraint the attack writes 6 meters (angle) or 7 (|V|) and
nothing changes. With it, both tests catch every draw from 1 deg or
0.02 pu. Below that it gets through: 0.25 deg moves the bus 7 injection by
about 3 MW, going by the 11.46 MW per degree from 11-25.

The J rise is much smaller than the 0.115 pu mismatch would suggest at a
0.001 pu std (115 sd). The estimator does not keep bus 8 where the
attacker put it and take the hit on the virtual rows; it moves bus 7 by
0.47 deg and bus 9 by 0.22 deg (noise-free run, 1 deg attack) to get S_7
back near zero, and the residual ends up mostly on the 4-7 and 7-9
transformer flows, then the bus 8 P meter and the 7-8 transformer flow.
Bus 8 itself goes to 1.21 deg, further than the attack. The virtual P row still has the largest r_N. With
the PMUs in z, buses 7 and 9 could not drift like that, so I think these
rates are a lower bound for the {2, 6, 9} placement, but I have not run
that.

The 0.001 pu virtual std is pandapower's constant, not something chosen
here.
