# grid-baddata-detect

Bad data and false data injection detection for power system state
estimation, on the IEEE 14-bus test case.

The project is simulation based: measurements come from a pandapower power
flow on the 14-bus case with Gaussian noise added. So far it runs
pandapower's WLS estimator on that measurement set and applies two bad
data tests: the chi-square test on J(x) and the largest normalised
residual test. To see both react to one gross error:

```
pip install -e .[dev]
python experiments/gross_error.py
```

`python experiments/detection_rate.py` repeats this over 200 noise draws
for a range of error sizes and plots how often each test fires, with both
thresholds set to 5% false alarms on a separate set of clean draws
(figures/detection_rate.png, about 17 minutes).

`python experiments/stealthy_attack.py` shifts the estimated angle of bus
14 through 10 meters without either test noticing, and compares that with
the linear a = Hc version and a random error of the same size
(figures/stealthy_attack.png, about 4 minutes).

`python experiments/pmu_check.py` adds one PMU voltage angle, out of the
attacker's reach, and shows it catches that attack only when it sits on
bus 14 itself (figures/pmu_check.png, about 3 minutes).

`python experiments/pmu_placement.py` finds how few protected PMUs leave
no state an attacker can move through the SCADA meters alone. With branch
currents it takes four, for example on buses 2, 6, 7 and 9, or three
(buses 2, 6 and 9) once the estimator knows that bus 7 has zero injection
(figures/pmu_placement.png, about 10 seconds).

`python experiments/zero_injection_attack.py` runs the one attack those
three PMUs leave open, on the angle and |V| of bus 8, through the
estimator with and without the bus 7 zero injection. Without it the tests
fire no more often than on clean data; with it both catch every draw from
a 1 degree or 0.02 pu shift (figures/zero_injection_attack.png, about 3
minutes).

HOW_IT_WORKS.md explains the method, what the simulation assumes and
where that flatters the results. NOTES.md has the dated log behind every
number.
