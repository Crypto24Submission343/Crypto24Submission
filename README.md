# Code for Heuristics in Appendix

Follow `README.md` in 

- AppendixB_Heu3
- AppendixE_Heu2

# Estimates for Kyber's QPE(W) in Section 3.1

Follow `README.md` in

- LWEEstimator

# Other Code

## Setup without Docker 

Follow the README files in

- Cost_Estimation
- fplll_experiments: Figures in appendices C, D, I
- AppendixE_Heu2: Code for generating cylinder pruning parameters for Kyber can also be found here.

## Setup with Docker

Run `sudo sh resetAndBuildDocker.sh`
Run `sudo sh runMyDocker.sh`

... or the few commands written in those files.

Dependencies required to run the experiments will be installed inside the container.
The experiments will then be run. The user will be left in a shell in the container. 
However, the results will be saved to `./costResults` and `./latticeExperiments` in the host machine.

Be aware: Running all of the code requires multiple hours, depending on the machine.
