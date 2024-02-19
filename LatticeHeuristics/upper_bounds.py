import subprocess
import collections
import json
import pathlib
from sage.all import log, factorial, is_even, ceil, sqrt
from functools import cache as memoize
from typing import List
from math import isnan
from lattices import gaussian_heuristic, vol_ball, gsa_basis_covol, RRR


import pickle
class MemoizeMutable:
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}
    def __call__(self, *args, **kwds):
        str = pickle.dumps(args, 1) + pickle.dumps(kwds, 1)
        if str not in self.memo: 
            # print("miss")  # DEBUG INFO
            self.memo[str] = self.fn(*args, **kwds)
        else:
            # print("hit")  # DEBUG INFO
            pass

        return self.memo[str]



def get_pruning_bound(dim, prec, target_pr):
    """
    Usage: ./test_pruner [options]
    List of options:
      -d <dimension>          (default=50)  Enumeration dimension
      -p <precision>          (default=212) Floating point precision
      -t <target_probability> (default=.51) Single enumeration target success probability
      -j                                    Get output in JSON format
    """
    script_path = pathlib.Path(__file__).parent.absolute()
    exec_path = pathlib.Path(f"{script_path}/../SubtreeDescendants/fplll/tests/").resolve()
    shellcmd = f'bash -c "cd {exec_path}; make test_pruner > /dev/null 2>&1; ./test_pruner -d {dim} -p {prec} -t {target_pr} -j; exit $?"'
    proc = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    out = collections.OrderedDict(json.loads(proc.stdout.decode('ascii').strip()))
    return out


def my_poly_integrate(P: List[RRR], len_P: int):
    for i in range(len_P-1, 0, -1):
        P[i] = P[i-1]/i
    P[0] = RRR(0.)
    return


def my_poly_eval(x: RRR, P: List[RRR], len_P: int):
    y = RRR(0.)
    pow = RRR(1.)
    for i in range(len_P):
        # y += P[i] * x^(i-1)
        y += P[i] * pow
        pow *= x
    return y


def my_poly_neg(P: List[RRR], len_P: int):
    for i in range(len_P):
        P[i] = -P[i]
    return


def VolumePolytope(Ri: List[float], P: List[RRR], access_even_radii_only=True):
    n = len(Ri)
    for i in range(1, len(P)):
        P[i] = RRR(0)
    P[0] = RRR(1)

    if access_even_radii_only:
        assert(is_even(n))
        step = -2
    else:
        step = -1

    for i in range(n-1, -1, step):
        ratio = RRR(Ri[i])/Ri[n-1]
        my_poly_integrate(P, n+1)
        newc = my_poly_eval(ratio**2, P, n+1)
        my_poly_neg(P, n+1)
        P[0] += newc

    rv = RRR(P[0])
    assert(rv > 0)
    return rv


@memoize
def mfactorial(x: int):
    return factorial(x)


def Pr_Chen13_level_cost(Ri, P: List[RRR], verbose=False, access_even_radii_only=True):
    n = len(Ri)
    if access_even_radii_only:
        assert(is_even(n))
        l = int(n/2)
    else:
        l = n
    fac = mfactorial(l)
    vol = RRR(VolumePolytope(Ri, P, access_even_radii_only))
    prod = fac * vol
    if verbose or prod <= 0:
        print("fac: ", fac, " (before Float64 cast)")
        print("vol: ", vol, " (before Float64 cast)")
        print("prod: ", prod, " (before Float64 cast)")
    rv = float(prod)
    assert(rv > 0)
    return rv


def eval_poly(ld: int, p: List[RRR], x: float) -> RRR:
    acc = RRR(0.0)
    xx = RRR(x)
    for i in range(ld, -1, -1):
        acc = acc * xx
        acc = acc + p[i]
    return acc

def integrate_poly(ld: int, p: List[RRR]):
    for i in range(ld, -1, -1):
        p[i + 1] = p[i] / RRR(i + 1)
    p[0] = RRR(0.)


def relative_volume(rd: int, b: List[RRR], factorial_tab=False):
    P: List[RRR] = [0. for _ in range(rd+2)]
    P[0]   = RRR(1)
    # ld = rd
    ld = 0
    for i in range(rd-1, -1, -1):
        integrate_poly(ld, P)
        ld += 1
        P[0] = RRR(-1.0) * eval_poly(ld, P, b[i] / b[rd - 1])

    fac = RRR(factorial(rd))
    res = RRR(P[0]) * fac

    if rd % 2 == 0:
        return res
    else:
        return -res


def Pr_Chen13_succ_prob_exact(Ri, P: List[RRR] = None, verbose: bool = False, fplll=False):
    """
    Algorithm 12 of Chen13.
    Assumes the enumeration radius = the norm of the target vector.
    It generalises the original analysis in GNR10 to tohe case where the dimension is odd.
    """
    n = len(Ri)

    _Ri = [Ri[i] / Ri[n-1] for i in range(n)]

    ceil_n_halfs = ceil(n/2)
    rp  = [_Ri[i-1] for i in range(1, 2 * ceil_n_halfs - 3 + 1, 2)]
    rpp = [_Ri[i-1] for i in range(2, 2 * ceil_n_halfs - 2 + 1, 2)]
    assert(len(rp) == ceil_n_halfs - 1)
    assert(len(rpp) == ceil_n_halfs - 1)
    if fplll:
        p_low  = relative_volume(len(rp), [r**2 for r in rp], factorial_tab=True)
        p_high = relative_volume(len(rpp), [r**2 for r in rpp], factorial_tab=True)
    else:
        fac = mfactorial(ceil_n_halfs - 1)
        if not P:
            P: List[RRR] = [0. for _ in range(n+2)]
        p_low  = fac * VolumePolytope( rp, P, access_even_radii_only=False)
        p_high = fac * VolumePolytope(rpp, P, access_even_radii_only=False)
    p_avg = (p_low + p_high)/2.
    return p_avg


def Pr_Chen13_succ_prob_exact_general(norm_v, Ri, P: List[RRR] = None):

    n = len(Ri)

    if not P:
        P: List[RRR] = [0. for _ in range(n+2)]

    if norm_v > Ri[-1]:
        return 0

    if norm_v == Ri[-1]:
        return Pr_Chen13_succ_prob_exact(Ri, P)

    # norm_v < Ri[-1]
    trim_Ri_norm_v = [min(norm_v, Ri[i]) for i in range(n)]
    return Pr_Chen13_succ_prob_exact(trim_Ri_norm_v, P)

def t(n, R, vol_L = 1.):
    v_n = vol_ball(n, 1)
    return RRR(R**n * v_n / vol_L / 2)

def kbar(n, R, R_p_D):
    return RRR((n/(n+1)) * (R_p_D*t(n, R_p_D) - R*t(n, R))/(t(n, R_p_D) - t(n, R)))

def Pr_Chen13_succ_prob_le(R, Ri):
    """
    Algorithm 13
    """
    n = len(Ri)
    p = 0
    u = 0
    s = [0 for _ in range(1001)]
    for i in range(1, 1001):
        s[i] = (i/1000)**(4/n) * R
        q = (1-p) * (1 - RRR(Pr_Chen13_succ_prob_exact_general(kbar(n, s[i-1], s[i]), Ri)
                                                            **(t(n, s[i])-t(n, s[i-1])
        )))
        p = p + q
        u = u + q * kbar(n, s[i-1], s[i])
    exp_R = u/p
    return float(exp_R), float(p)

# def gsa_basis_covol(dim: int, k: int, V: float = 1.):
#     # this is for a lattice of volume 1
#     alpha = gsa_alpha(dim)
#     prod = 1.
#     b1norm = alpha ** (-(dim-1)/2.)
#     for j in range(1, k+1):
#         i = dim-j+1
#         bistar = alpha ** (i-1) * b1norm * V**(1/dim)
#         prod *= bistar
#     assert(prod > 0)
#     return prod


@MemoizeMutable
def CylinderPruningUB_normalise_Ri(n: int, Ri: List[RRR], V: float = 1., eps: float = 10**-10) -> List[RRR]:
    # normalise pruning radii
    assert(n == len(Ri))
    R = gaussian_heuristic(V, n)
    _Ri = [R * Ri[i] / Ri[n-1] for i in range(n)]
    if abs(_Ri[n-1] - R) < eps:
        return _Ri
    else:
        print("Error in normalisation, abs(Ri[n-1] - R) =", abs(_Ri[n-1] - R))
    assert(abs(_Ri[n-1] - R) < eps)


@MemoizeMutable
def CylinderPruningUB_log_Hk(k: int, n: int, Ri: List[float], P: List[RRR] = None, normalise_Ri=True, verbose=False) -> RRR:
    assert(len(Ri) == n)
    assert(k > 0)
    assert(k <= n)

    two = RRR(2.)

    if not P:
        P: List[RRR] = [0. for _ in range(n+2)]

    if normalise_Ri:
        # normalise pruning radii
        _Ri = CylinderPruningUB_normalise_Ri(n, Ri)
    else:
        _Ri = Ri

    # we are assuming unit volume
    covol = gsa_basis_covol(n, k)
    if k == 1:
        # the 1-ball of radius _Ri[0] is the [-Ri[0], Ri[0]] interval
        vol_k_ball_R1 = 2 * _Ri[0]
        return log(vol_k_ball_R1, 2) - log(covol, 2)

    if is_even(k):
        vol_k_ball = vol_ball(k, _Ri[k-1])
        even_Ri = [Ri[i] for i in range(1, k, 2)]
        odd_Ri  = [Ri[i] for i in range(0, k, 2)]
        upper_estimate_pr = Pr_Chen13_level_cost(even_Ri, P, access_even_radii_only=False)
        lower_estimate_pr = Pr_Chen13_level_cost(odd_Ri,  P, access_even_radii_only=False)
        avg_pr = (upper_estimate_pr + lower_estimate_pr)/two
        return log(vol_k_ball, 2) + log(avg_pr, 2) - log(covol, 2)
        # return log(vol_k_ball, 2) + log(avg_pr, 2)
        # return log(vol_k_ball, 2) - log(covol, 2)
        # return log(avg_pr, 2) - log(covol, 2)
        # return log(vol_k_ball, 2)
        # return log(avg_pr, 2)

    if k == n:
        assert(not is_even(k))
        # we are doing extreme pruning, the last level likely is empty
        # log2(0) ≈ -9999
        # about 1/pr points in the last level
        return log(Pr_Chen13_succ_prob_exact(Ri), 2)

    # k != 1, k is odd, k != n
    # average level above and below
    # Note: either called with Ri already normalised, or Ri was normalised at the beginning of this function
    above = CylinderPruningUB_log_Hk(k-1, n, Ri)
    below = CylinderPruningUB_log_Hk(k+1, n, Ri)
    if verbose:
        print("above", float(above))
        print("below", float(below))
    # avg = (two ** above + two ** below)/two
    avg = two**((above+below)/2) # same interpolation as fplll, smoother curve
    return log(avg, 2)


def CylinderPruningUB_log_cost(n, Ri):
    two = RRR(2.)
    P: List[RRR] = [0. for _ in range(n+2)]

    # normalise pruning radii
    _Ri = CylinderPruningUB_normalise_Ri(n, Ri)

    N = RRR(0.)
    for k in range(1, n+1):
        N += two ** CylinderPruningUB_log_Hk(k, n, Ri)

    pr = Pr_Chen13_succ_prob_exact(Ri, P)

    single_enum_log_cost = float(log(N, 2)) - 1. # -1 due to SVP symmetry
    total_enum_log_cost = single_enum_log_cost - float(log(pr, 2))
    return single_enum_log_cost, total_enum_log_cost


##

