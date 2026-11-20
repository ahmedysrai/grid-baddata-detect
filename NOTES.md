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
