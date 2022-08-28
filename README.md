# 

## The MWD class


## The Trigger class


---

### Comparing against `MWD.m`
This code is meant to give the same result as `MWD.m` by James Lawson. To check this, when running `MWD.py` as main, a set of local variables can be written to a file, which can then be compared with a reference file written by `MWD.m`. 

To generate the reference file, add the following lines to the end of MWD.m:
```(language=matlab)
tfilt_padded = [zeros(1,cfd_delay) tfilt'];
outfileID = fopen("./original_variables.tsv", 'w');
fprintf(outfileID,'%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%3d\t%3d\n',...
        [D; MA; MWD; T; tfilt_padded; BL; trigger_event*111; sample_time*111]);
```

The function `_write_local_variables` in `MWD.py` generates a similar file. 

The function `_test_local_variables` in `MWD.py` compares the two files and prints the largest discrepancy for each variable.
**Note that that the input parameters are not compared, meaning the user has to make sure they are the same.**

For convenience, reference files are provided for the test traces in [](test_traces).
