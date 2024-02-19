# Appendix E experiments

This code is based on the code from fplll_experiments.
It does however generate significantly larger stat files, in the GB for n = 60.
Hence we separate it to avoid unnecessary computation.

Follow the same preparation instructions as for fplll_experiments.
Then run
```
make run_jensen_experiments
```
to reproduce the results from scratch, or
```
bash -c "cd pruned; tar xf stats.tar.gz"
make plot_jensen_experiments
```
to regenerate the plots from the preexisting data.

The stat files and plots presented in the paper are already stored.

## Results

The plots generated and used for the paper are in `pruned-paper`.
We notice that due to a typo in our code, the plots in the paper are suboptimal.
The correct plots, with results better matching our predictions, can be found in `pruned`.

# Pruning parameters

The `batch_pruning_radii_gen.sh` file can be used to generate cylinder pruning parameters for Kyber-relevant dimensions.
Before running, `test_pruner` needs to be compiled:

```bash
cd fplll/tests
make test_pruner
```

Running `batch_pruning_radii_gen.sh` will run all processes in parallel. They can take various hours.
