The code in this directory is used to replicate our enumeration tree size estimates for Kyber in Section 3.1.

1. install docker, and run:
sudo docker pull sagemath/sagemath

2. from within the directory containing this file, run:
sudo docker run -v $(pwd)/..:/mnt -it sagemath/sagemath bash

3. build aflk library:
cd /mnt/LatticeHeuristics/alfk
python3 build_alfk.py

4. run estimates
sage /mnt/LWEEstimator/our_kyber_estimates.py