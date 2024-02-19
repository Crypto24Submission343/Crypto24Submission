from sage.all import cached_function, RealField, RR, ZZ, pi, e, gamma, sqrt

RRR = RealField(1000)


@cached_function
def gaussian_heuristic(v, n):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """ 
    from math import pi, sqrt, e
    v = RRR(v)
    n = ZZ(n)
    pi_part = ((pi * n)**(1/(2*n)))
    sqrt_part = sqrt(n/(2 * pi * e))
    v_part = (v ** (1/n))
    # print(f"v {v}, pi_part {pi_part}, sqrt_part {sqrt_part}, v_part {v_part}")
    return pi_part * sqrt_part * v_part


@cached_function
def stirling_gamma(x):
    x = RRR(x)
    if x < 10000:
        return gamma(x)
    else:
        raise ValueError("too large")
        _gamma = ((2*pi*(x-1))**RRR(.5)) * (((x-1)/e)**(x-1))
        return _gamma


@cached_function
def vol_ball(dim, rad):
    dim, rad = ZZ(dim), RRR(rad)
    rad_part = rad**dim
    pi_part = pi**(dim/2)
    gamma_part = stirling_gamma(dim/2 + 1)
    # print(f"vol_ball rad {rad_part}, pi {pi_part}, gamma {gamma_part}")
    vol = rad_part * pi_part / gamma_part
    return vol


@cached_function
def deltaf(beta):
    """
    Taken from https://github.com/malb/lattice-estimator/blob/f38155c863c8a1fb8da83a10664685d3496fe219/estimator/reduction.py#L12

    Compute δ from block size β without enforcing β ∈ ZZ.
    δ for β ≤ 40 were computed as follows:
    ```
    # -*- coding: utf-8 -*-
    from fpylll import BKZ, IntegerMatrix
    from multiprocessing import Pool
    from sage.all import mean, sqrt, exp, log, cputime
    d, trials = 320, 32
    def f((A, beta)):
        par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, flags=BKZ.AUTO_ABORT)
        q = A[-1, -1]
        d = A.nrows
        t = cputime()
        A = BKZ.reduction(A, par, float_type="dd")
        t = cputime(t)
        return t, exp(log(A[0].norm()/sqrt(q).n())/d)
    if __name__ == '__main__':
        for beta in (5, 10, 15, 20, 25, 28, 30, 35, 40):
            delta = []
            t = []
            i = 0
                while i < trials:
                threads = int(open("delta.nthreads").read()) # make sure this file exists
                pool = Pool(threads)
                A = [(IntegerMatrix.random(d, "qary", beta=d//2, bits=50), beta) for j in range(threads)]
                for (t_, delta_) in pool.imap_unordered(f, A):
                    t.append(t_)
                    delta.append(delta_)
                i += threads
                print u"β: %2d, δ_0: %.5f, time: %5.1fs, (%2d,%2d)"%(beta, mean(delta), mean(t), i, threads)
            print
    ```
    """
    small = (
        (2, 1.02190),  # noqa
        (5, 1.01862),  # noqa
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )

    if beta <= 2:
        return RR(1.0219)
    elif beta < 40:
        for i in range(1, len(small)):
            if small[i][0] > beta:
                return RR(small[i - 1][1])
    elif beta == 40:
        return RR(small[-1][1])
    else:
        return RR(beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


@cached_function
def gsa_alpha(beta):
    delta = RRR(deltaf(beta))
    return delta ** -2   


@cached_function
def gsa_basis_covol(dim: int, k: int, V: float = 1.):
    # this is for a lattice of volume 1
    alpha = gsa_alpha(dim)
    prod = 1.
    b1norm = alpha ** (-(dim-1)/2.)
    for j in range(1, k+1):
        i = dim-j+1
        bistar = alpha ** (i-1) * b1norm * V**(1/dim)
        prod *= bistar
    assert(prod > 0)
    return prod
