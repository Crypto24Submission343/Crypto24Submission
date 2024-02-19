Instructions to run this code.
==============================

1. install docker, and run:
sudo docker pull sagemath/sagemath

2. from within the directory containing this file, run:
sudo docker run -v $(pwd):/mnt -it sagemath/sagemath bash 

3. inside the bash session that opens run once:
pip install tqdm seaborn

# Data collection

Inside the Docker instance run the following.

# No pruning experiments:
```bash
sudo bash -c "sage /mnt/covariance.py -m 46 -M 50 -gh 1.05 -t 30 -tpe 2  >> /mnt/cov_res.txt"
sudo bash -c "sage /mnt/covariance.py -m 52 -M 54 -gh 1.05 -t 30 -tpe 2 --prec 60 >> /mnt/cov_res.txt"
```

# Linear pruning experiments:
NOTE: requires changing a couple of lines in covariance.py (search "For linear pruning experiments")
```bash
sudo bash -c "sage /mnt/covariance.py -m 40 -M 66 -gh 1.05 -t 30 -tpe 2 --bits 5 >> /mnt/lin_pru.txt"
```

Note that for running on smaller dimensions, a more relaxed gaussian heuristic factor than 1.05 may be required (say, 1.15 or 1.20) to find short vectors.
Keeping such a large factor throughout will significantly slow down experiments at larger dimension.

The output will be stored in `cov_res.txt` in the directory where this file is contained.

# Plot generation

Inside the Docker instance run:

```bash
sage /mnt/plot_cov_res.py -i /mnt/cov_res_np_46_50.txt -o /mnt/covariance/np
sage /mnt/plot_cov_res.py -i /mnt/cov_res_np_52_54.txt -o /mnt/covariance/np
sage /mnt/plot_cov_res.py -i /mnt/lin_pru.txt -o /mnt/covariance/lp
```
