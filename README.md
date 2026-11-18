# grid-baddata-detect

Bad data and false data injection detection for power system state
estimation, on the IEEE 14-bus test case.

The project is simulation based: measurements come from a pandapower power
flow on the 14-bus case with Gaussian noise added. So far it runs
pandapower's WLS estimator on that measurement set and computes the
chi-square statistic J(x). To see J react to one gross error:

```
pip install -e .[dev]
python experiments/gross_error.py
```

NOTES.md has the results so far.
