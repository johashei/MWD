# 

## The MWD class
This class implements the moving window deconvolution. Assuming a signal 
```math
P(t) = A*\mathrm{e}^{-t/\tau} \quad \mathrm{for } t\geq 0, 
```
the pulse heigh $$A$$ can be calculated by
```math
A = P(t) + \frac{1}{\tau} \int_{t'=0}^{t} P(t')\mathrm{d}t' + P(t).
```
This expression is discretised as
```math
A_i = P_i - P_{i-M} + \frac{1}{\tau}\sum_{k=i-M}^{i-1}P_k
```

## The Trigger class


---

### Comparing against `MWD.m`
This code is meant to give the same result as `MWD.m` by James Lawson. To check this, when running `MWD.py` as main, a set of local variables can be written to a file, which can then be compared with a reference file written by `MWD.m`. 

To generate the reference file, add the following lines to the end of MWD.m:
```matlab
tfilt_padded = [zeros(1,cfd_delay) tfilt'];
outfileID = fopen("./matlab_local_variables.tsv", 'w');
fprintf(outfileID,'%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%3d\t%3d\n',...
        [D; MA; MWD; T; tfilt_padded; BL; trigger_event*111; sample_time*111]);
```

The function `_write_local_variables` in `MWD.py` generates a similar file. 

The function `_test_local_variables` in `MWD.py` compares the two files and prints the largest discrepancy for each variable.
**Note that that the input parameters are not compared, meaning the user has to make sure they are the same.**

For convenience, reference files are provided for the test traces in [test_traces](./test_traces).
These were generated using the following parameters:
```python
## CFD parameters
cfd_delay = 8
cfd_threshold = 150
glitch_filter_threshold = 75
## MWD parameters 
trapezoid_length = 600
rise_time = 450
decay_time = 5e3
## Baseline estimation parameters
baseline_length = trapezoid_length + rise_time + 110
baseline_fit_window = 100 
```
