"""
Instructions to run this code.
==============================

1. install docker, and run:
sudo docker pull sagemath/sagemath

2. from within the directory containing this file, run:
sudo docker run -v $(pwd):/mnt -it sagemath/sagemath bash 

3. inside the bash session that opens run once:
pip install tqdm

4. finally run:
sudo bash -c "sage /mnt/covariance.py -m 44 -M 50 -gh 1.05 -t 30 -tpe 2  >> /mnt/cov_res.txt"

5. For larger experiments run with higher-than-double precision (slower)
sudo bash -c "sage /mnt/covariance.py -m 52 -M 52 -gh 1.05 -t 30 -tpe 2 --prec 60 >> /mnt/cov_res.txt"

Note that for running on smaller dimensions, a more relaxed gaussian heuristic factor than 1.05 may be required (say, 1.15 or 1.20) to find short vectors.
Keeping such a large factor throughout will significantly slow down experiments at larger dimension.

The output will be stored in `cov_res.txt` in the directory where this file is contained.
"""


from fpylll import IntegerMatrix, Enumeration, EnumerationError, GSO, BKZ, load_strategies_json, LLL, FPLLL
import fpylll
from sage.all import matrix, sqrt, log, randint
import numpy
from numpy import mean
import math


# we define variance() due to deprecation by sage
def variance(l):
    return numpy.var(l, ddof=1) # ddof=1 for unbiased estimator

def covariance(l):
    """
        l = list of samples of n random variables, where n = len(l[i]) for all i < len(l)
        outputs a covariance matrix where C[i,j] = Cov(Xi, Xj)
    """
    x = numpy.array(l)
    return numpy.cov(x, rowvar=False)

def gaussian_heuristic(v, n, verbose=True):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """
    return ((math.pi * n)**(1./(2.*n))) * math.sqrt(n/(2 * math.pi * math.e)) * (v ** (1./n))


from sage.all import RR
def vol_n_ball(n, R):
    try:
        return float(R**n * math.pi**(n/2) / math.gamma(n/2+1))
    except:
        return float(RR(R)**n * RR(math.pi)**(n/2) / RR(math.gamma(n/2+1)))


# DEF_GH_FACTOR = 1.10
sampled_primes = {}


from sage.all import random_prime, vector, identity_matrix, randrange
def random_cn11_lattice(n, q):
    m = identity_matrix(n-1).augment(vector([randrange(q) for _ in range(n-1)])).stack(vector([0] * (n-1) + [q]))
    B = IntegerMatrix.from_matrix(m)
    return B


def experiment(n_and_prng_seed_and_gh_factor_and_threads_and_prec_and_bit_factor, max_sols=2**31-1):
    """
    This function generates a random basis for some kind of lattice, and runs enumeration over it.
    It then reports the number of solutions found, and the number of nodes explored as part of enumeration.

    It should be noticed that if no solutions are found, it rases an exception. If more solutions than max_sols is found,
    fpylll will return at most max_sols solutions. This means that for lattices of high enough dimension, one ends up getting
    all the calls return exactly max_sols vectors, and hence the data is invalid. For this reason, in such case an exception is raised.
    """
    n, prng_seed, gh_factor, threads, prec, bit_factor = n_and_prng_seed_and_gh_factor_and_threads_and_prec_and_bit_factor
    prng_seed = int.from_bytes(prng_seed, byteorder="big")

    solved = False
    cnt = 0
    while not solved: # at small dimension, no points may be in the enum radius. In that case, retry with a new prg seed
        cnt += 1
        if cnt > 20:
            raise Exception('inf loop')

        B = IntegerMatrix(n, n)

        # generate a random basis using fpylll
        fpylll.util.set_random_seed(prng_seed)
        bits  = bit_factor*n
        q = sampled_primes[bits]
        # q = large_enough_prime
        with seed(prng_seed):
            B = random_cn11_lattice(n, q)
            pass

        # generate a random basis using sagemath
        from sage.crypto import lattice
        from sage.all import next_prime
        with seed(prng_seed):
            # B = IntegerMatrix.from_matrix(lattice.gen_lattice("random", 1, n, next_prime(500)))
            pass

        if prec <= 53:
            GS = GSO.Mat(B)
        else:
            FPLLL.set_precision(prec)
            GS = GSO.Mat(B, float_type="mpfr")
        GS.update_gso()
        L = LLL.Reduction(GS)
        L()

        M = matrix(n)
        B.to_matrix(M)
        covol = abs(M.determinant())

        # different bases will have different covolumes, rather than unit covolume. to address this, pick the radius of the ball so that
        # vol(ball) = covol. Below, we do exactly this, except we also scale up the ball so that it's a bit larger than the covolume (increasing
        # the chance that we don't get balls that happen to be too small for any vectors to be in it for some unlucky basis)

        gh = gaussian_heuristic(covol, n)
        R_sqrd = (gh_factor * gh)**2
        vol_ball = vol_n_ball(n, sqrt(R_sqrd))
        try:
            FPLLL.set_threads(threads)
            E = Enumeration(GS, nr_solutions=max_sols)
            try:
                # For linear pruning experiments
                # p = fpylll.Pruning.LinearPruningParams(n, n)
                # enumres = E.enumerate(0, n, R_sqrd, 0, pruning=p.coefficients)

                # For no-priuning experiments
                enumres = E.enumerate(0, n, R_sqrd, 0)

                nsols = len(enumres)
                if nsols == max_sols:
                    raise(Exception("more than max_sols found"))
                nodes  = E.get_nodes()
                H_n_minus_k = [E.get_nodes(level=n_minus_k) for n_minus_k in range(n)] # I think H_n_minus_k(0) corresponds to k = n and H_n_minus_k(n-1) corresponds to k = 1
                H_n_minus_k_inv = [1/_ for _ in H_n_minus_k]
                # import pdb; pdb.set_trace()
                svpnorm = sqrt(min([sqr_norm for sqr_norm, _ in enumres]))
            except EnumerationError as e:
                raise(e)
                nsols = 0
                nodes = 0
                svpnorm = -1
            solved = True
        except:
            pass
        prng_seed = randint(0, 1000000000)

    return { 'nsols': nsols, 'vol': vol_ball/covol, 'nodes': nodes, 'H_n_minus_k_with_inv': H_n_minus_k + H_n_minus_k_inv, 'svnorm': svpnorm, 'covol': covol }


import tqdm
import multiprocessing
nproc = int(0.8 * multiprocessing.cpu_count()) - 2
print(f"Using {nproc} cpus")


def main(tries=4, min_n=40, max_n=50, gh_factor=1.1, threads_per_enum=1, prec=0, bit_factor=10):
    nsols = {}
    n = min_n
    try:
        # tab_format = '| {:2} | {:>18} | {:>18} | {:>18} | {:>18} | {:>18} | {:>18} | {:>18} | {:>18} | {:>18} |{:>18} |'
        # print(tab_format.format("n", "exp(nsols)/λ", "stddev(nsols)/λ", "var(nsols)/λ", "raw2(nsols)/raw2(poi)", "raw3(nsols)/raw3(poi)", "exp nodes", "stddev nodes", "var nodes", "λ", f"{DEF_GH_FACTOR}^n"))
        while n <= max_n:
            nsols[n] = {'nsols': [], 'nodes': [], 'H_n_minus_k_with_inv': [], 'vol': -1, 'svnorm': -1}

            # for _ in tqdm.tqdm(range(tries), position=0, desc=f"n: {n}", leave=False):
            #     rv = experiment(n)
            #     nsols[n]['nsols'].append(rv['nsols'])
            #     nsols[n]['vol'] = rv['vol']

            with multiprocessing.Pool(nproc) as pool:
                for rv in tqdm.tqdm(pool.imap_unordered(experiment, list((n, bytes(f"{n}-{prng_seed}", "ascii"), gh_factor, threads_per_enum, prec, bit_factor) for prng_seed in range(1, tries+1))), total=tries, leave=False):
                    nsols[n]['H_n_minus_k_with_inv'].append(rv['H_n_minus_k_with_inv'])
                    nsols[n]['nsols'].append(rv['nsols'])
                    nsols[n]['nodes'].append(rv['nodes'])
                    nsols[n]['vol'] = rv['vol']
                pool.close()
                pool.join()

            print("\n\n\n\n\n\n\n\n")
            print("=============================")
            print(f"tries={tries}, n={n}, gh_factor={gh_factor}, threads_per_enum={threads_per_enum}, prec={prec}, bit_factor={bit_factor}")
            print("\n\n\n\n\n\n\n\n")

            print(n, len(nsols[n]['H_n_minus_k_with_inv']), nsols[n]['nsols'])
            C = covariance(nsols[n]['H_n_minus_k_with_inv'])

            # # gets variances of Hk
            # for i in range(n):
            #     # print("n-k =", i)
            #     print("k =", n-i)
            #     print(C[i,i])
            #     print(variance([v[i] for v in nsols[n]['H_n_minus_k_with_inv']]))
            #     print()

            means = [
                mean([v[i] for v in nsols[n]['H_n_minus_k_with_inv']])
                for i in range(2*n)
            ]

            for k in range(1, n+1):
                if k == int(n/2):
                    print("\n\n\n")
                for h in range(1, n+1-k):
                    print("k: %02d" % k, end="; ")
                    print("k+h: %02d" % (k+h), end="; ")
                    n_k = n-k
                    n_k_h = n-k-h
                    print("Cov(1/|Z_k|, |Z_{k+h}|):", end="\t")
                    print("%0.20f" % C[n_k+n, n_k_h], end=";\t")

                    print("100 * |Cov(1/|Z_k|, |Z_{k+h}|)| / (H_{k+h}/H_k):", end="\t")
                    H_k = means[n_k]
                    H_kh = means[n_k_h]
                    print("%0.20f" % (100 * abs(C[n_k+n, n_k_h]) / (H_kh/H_k)), end=";\t")

                    Exp_ratio = C[n_k+n, n_k_h] + means[n_k_h] * means[n_k+n]
                    # print("E[X/Y]:", end="\t")
                    # print("%0.6f" % (Exp_ratio), end="\t")

                    print("(E[X]/E[Y]) / E[X/Y]:", end="\t")
                    print("%0.20f" % ((H_kh/means[n_k]) / Exp_ratio), end=";\t")
                    #
                    print("(E[X]=E[H_kh]:", end="\t")
                    print("%0.20f" % (H_kh), end=";\t")
                    #
                    print("E[Y]=E[H_k]:", end="\t")
                    print("%0.20f" % (H_k), end=";\t")
                    #
                    print("E[X/Y]=E[H_kh/H_k]:", end="\t")
                    print("%0.20f" % (Exp_ratio), end=";\t")
                    print()
                    # print("%0.6f" % ((H_kh*means[n_k+n]) / Exp_ratio)) ## roughe edit, ignores Jensen

            # poisson_param = (nsols[n]['vol']/2)
            # ghf = DEF_GH_FACTOR
            # print(tab_format.format(
            #     n,
            #     "%.2f" % (mean(nsols[n]['nsols'])/poisson_param),
            #     "%.1f" % sqrt(variance(nsols[n]['nsols'])/poisson_param),
            #     "%.2f" % (variance(nsols[n]['nsols'])/poisson_param),
            #     "%.2f" % (raw_k_moment_estimator(2, nsols[n]['nsols'])/poisson_k_moment(2, poisson_param)),
            #     "%.2f" % (raw_k_moment_estimator(3, nsols[n]['nsols'])/poisson_k_moment(3, poisson_param)),
            #     "%.1f" % mean(nsols[n]['nodes']),
            #     "%.1f" % sqrt(variance(nsols[n]['nodes'])),
            #     "%.1f" % variance(nsols[n]['nodes']),
            #     "%.1f" % poisson_param, # poisson parameter
            #     "%.2f" % (ghf**n)
            # ))
            n += 2
    except KeyboardInterrupt:
        pass

    return nsols


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-n', type=int, default=40)
    parser.add_argument('-M', '--max-n', type=int, default=50)
    parser.add_argument('-t', '--tries', type=int, default=4)
    parser.add_argument('-tpe', '--threads_per_enum', type=int, default=2)
    parser.add_argument('-gh', type=float, default=1.1)
    parser.add_argument('--prec', type=int, default=0)
    parser.add_argument('--bits', type=int, default=10)
    args = parser.parse_args()

    min_n = args.min_n
    max_n = args.max_n
    gh_factor = args.gh
    tries = args.tries
    threads_per_enum = args.threads_per_enum
    prec = args.prec
    bit_factor = args.bits

    # compute primes for CN11 random lattice generation
    from sage.all import seed
    with seed(0xdeadbeef):
        for n in range(1, args.max_n+1):
            bits  = bit_factor*n
            sampled_primes[bits] = random_prime(2**bits)
        large_enough_prime = sampled_primes[bit_factor*args.max_n]

    # experiment((20, bytes(0xdeadbeef)))
    main(tries=tries, min_n=min_n, max_n=max_n, gh_factor=gh_factor, threads_per_enum=threads_per_enum, prec=prec, bit_factor=bit_factor)
