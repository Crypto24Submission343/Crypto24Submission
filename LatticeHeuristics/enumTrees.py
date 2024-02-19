from sage.all import save, sqrt, log, prod, RR, pi, gamma, e, RealField, cached_function, ZZ, floor, infinity
from sage.all import beta as beta_function
from bkz_simulators.sim import Sim, LLLProfile
from alfk import _alfk
from lattices import gaussian_heuristic, vol_ball, gsa_alpha, deltaf, gsa_basis_covol, RRR, stirling_gamma
import collections
import functools
from typing import List
import pathlib
import json
import math
from upper_bounds import get_pruning_bound, Pr_Chen13_succ_prob_exact, CylinderPruningUB_log_cost, CylinderPruningUB_log_Hk

LOG2_PI = log(pi, 2)

RRRR = RealField(10000)


def fit_curve(fun, verbose=False, incl_logn=True, given_nlogn=None):
    from sage.all import var, find_fit
    dim = []
    nodes = []
    prob = []
    # for _dim in [250, 500, 1000, 1500]: # TODO: return to 1001
    gap = 300 # for fine graine, 30
    while len(dim) < (4 if incl_logn else 3):
        dim = []
        nodes = []
        prob = []
        for _dim in range(100, 1001, gap): # TODO: return to 1001
            if verbose: print("doing dim", _dim)
            try:
                _nodes, _prob = fun(_dim)
                _nodes = RR(_nodes)
            except:
                _nodes = RR(fun(_dim))
                _prob = 1
            if not _nodes.is_NaN():
                dim.append(_dim)
                nodes.append(_nodes)
                prob.append(_prob)
        T = list(zip(dim, nodes)) 
        gap -= 20

    if incl_logn:
        # fit assuming logn in the exponent too
        a, b, c, d, beta = var("a, b, c, d, beta") 
        if given_nlogn:
            a = given_nlogn
        f = a * beta * log(beta, 2.0) + b * beta + c * log(beta, 2.0) + d
        f = f.function(beta) 
        out = find_fit(T, f, solution_dict=True)
        if given_nlogn:
            return { 'nl2': a, 'n': out[b], 'l2': out[c], 'c': out[d] }, prob
        else:
            return { 'nl2': out[a], 'n': out[b], 'l2': out[c], 'c': out[d] }, prob
    else:
        a, b, c, beta = var("a, b, c, beta") 
        if given_nlogn:
            a = given_nlogn
        f = a * beta * log(beta, 2.0) + b * beta + c 
        f = f.function(beta) 
        out = find_fit(T, f, solution_dict=True)
        if given_nlogn:
            return { 'nl2': a, 'n': out[b], 'l2': 0, 'c': out[c] }, prob
        else:
            return { 'nl2': out[a], 'n': out[b], 'l2': 0, 'c': out[c] }, prob


def plot_asymptotic(exponent, filename=None):
    from sage.all import line, save
    f = lambda x: x * log(x,2) * exponent['nl2'] + x * exponent['n'] + log(x,2) * exponent['l2'] + exponent['c']
    g = line([(x, f(x)) for x in range(50, 500)])
    if filename:
        save(g, filename)
    return g


class NoPruning:

    @staticmethod
    @cached_function
    def log_Hk(k, n, exact=True, R=None, V=None, **kwargs):
        """
            :param R: enumeration radius
            :param V: volume of lattice
        """
        k, n = ZZ(k), ZZ(n)

        if k == 0:
            # only one node on root level
            return 0

        alpha = gsa_alpha(n)
        if exact:
            if not V:
                # assume covol of the full lattice = 1
                V = 1
            if not R:
                R = gaussian_heuristic(V, n)
            # print("enumTrees > R = ", R)
            # print("enumTrees > GH = ", gaussian_heuristic(1, n))
            # print("enumTrees > log V = ", float(log(V,2)), "    V^(1/n) = ", V**(1/n))
            b1norm = alpha**(-(n-1)/2)
            vol = vol_ball(k, R)
            prod = 1
            for j in range(1, k+1):
                i = n-j+1
                bistar = alpha**(i-1) * b1norm
                if V != 1:
                    bistar = bistar * V**(1/n) # adjust volume if passed in input
                prod *= bistar
            return log(vol,2) - log(prod,2)
        else:
            if R:
                raise ValueError("inexact non-R computation of log_Hk is not implemented")
            exp_2 = k * log(n/k, 2) / 2 + (k/n - 1) * (LOG2_PI + log(n, 2)) / 2
            exp_a = k * (k - n) / 2
            return exp_2 + exp_a * log(alpha, 2)

    @staticmethod
    @cached_function
    def log_Hk_single_enum(k, n, exact=True, R=None, **kwargs):
        return NoPruning.log_Hk(k, n, exact=True, R=None, **kwargs)

    @staticmethod
    @cached_function
    def log_cost(n, **kwargs):
        n = ZZ(n)
        N = RRRR(0)
        for k in range(1, n+1):
            N += RRRR(2) ** NoPruning.log_Hk(k, n)
        return log(N, 2)-1

    @staticmethod
    @cached_function
    def log_single_enum(n, **kwargs):
        return NoPruning.log_cost(n)

    @staticmethod
    @cached_function
    def log_Skh(norm_g, k, h, n, **kwargs):
        """ Assumes unit volume.
        returns the number of nodes on the (k+h)'s enumeration tree level that descend directly from g
         """
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        norm_g = RR(norm_g)
        alpha = gsa_alpha(n)
        # assume unit volume lattice
        R = gaussian_heuristic(1, n)
        b1norm = alpha**(-(n-1)/2)
        vol = vol_ball(h, sqrt(R**2 - norm_g**2))

        prod = 1
        for j in range(1, h+1):
            i = n-k-j+1
            bistar = alpha**(i-1) * b1norm
            prod *= bistar
        return log(vol,2) - log(prod,2)

    @staticmethod
    @cached_function
    def log_subtree_N_kD(norm_g, k, D, n, **kwargs):
        """
        Returns the size of the subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D by a vector on level k of norm norm_g
        assumes unit olume
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        norm_g = RR(norm_g)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** NoPruning.log_Skh(norm_g, k, h, n)
        return log(N, 2)

    @staticmethod
    @cached_function
    def log_avg_Skh_ub(k, h, n, **kwargs):
        """ upper-bounds S_{k,h}(g) by replacing g_norm^2 with the average square norm for g on level k 
        Assumes unit volume
        returns the number of nodes on the (k+h)'s enumeration tree level that descend directly from g for average g on level k """
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        # assume unit volume lattice
        R = gaussian_heuristic(1, n)
        avg_norm_g = R * sqrt(k / (k+2))
        return NoPruning.log_Skh(avg_norm_g, k, h, n)

    @staticmethod
    @cached_function
    def avg_children_Ck(k, n, **kwargs):
        return 2**NoPruning.log_avg_Skh(k, 1, n, **kwargs)

    @staticmethod
    @cached_function
    def log_avg_subtree_N_kD_ub(k, D, n, **kwargs):
        """
        Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
        assumes unit olume
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** NoPruning.log_avg_Skh_ub(k, h, n)
        return log(N, 2)

    @staticmethod
    @cached_function
    def log_avg_Skh(k, h, n, new_analysis=True, **kwargs):
        """ returns the number of nodes on the (k+h)'s enumeration tree level that descend directly from g for average g on level k """
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        if new_analysis:
            return NoPruning.log_Hk(k+h, n) - NoPruning.log_Hk(k, n)
        else:
            # assume unit volume lattice
            alpha = gsa_alpha(n)
            b1norm = alpha**(-(n-1)/2)
            R = gaussian_heuristic(1, n)
            V_k_ph = vol_ball(k+h, R)
            V_k = vol_ball(k, R)
            prod = 1
            for j in range(1, h+1):
                i = n-k-j+1
                bistar = alpha**(i-1) * b1norm
                prod *= bistar
            return log(V_k_ph, 2) - log(V_k, 2) - log(prod, 2)

    @staticmethod
    @cached_function
    def log_avg_Skh_single_enum(k, h, n, new_analysis=True, **kwargs):
        return NoPruning.log_avg_Skh(k, h, n, new_analysis=new_analysis, **kwargs)

    @staticmethod
    @cached_function
    def log_avg_subtree_N_kD(k, D, n, **kwargs):
        """
        Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
        assumes unit volume
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** NoPruning.log_avg_Skh(k, h, n)
        return log(N, 2)


class LinearPruning:

    @staticmethod
    @cached_function
    def log_Hk(k, n, bound="lower", **kwargs):
        """
        Hk for a single tree. Would succeed with prob 1/n.
        """
        k, n = ZZ(k), ZZ(n)

        if k == 0:
            # only one node on root level
            return 0

        if bound == "upper":
            poly_factor = 1
        elif bound == "lower":
            poly_factor = 1/k
        else:
            raise ValueError("`bound` can only be \"lower\" or \"upper\"")

        exp_no_pruning = NoPruning.log_Hk(k, n)
        exp_lin_pruning_2 = k * log(k / n, 2) / 2 + log(poly_factor, 2)
        return exp_no_pruning + exp_lin_pruning_2

    @staticmethod
    @cached_function
    def log_Hk_single_enum(k, n, bound="lower", **kwargs):
        return LinearPruning.log_Hk(k, n, bound=bound, **kwargs)

    @staticmethod
    @cached_function
    def log_cost(n, bound="lower", **kwargs):
        """
        Total size of a single tree. Would succeed with prob 1/n.
        """
        n = ZZ(n)
        N = RRRR(0)
        for k in range(1, n+1):
            N += RRRR(2) ** LinearPruning.log_Hk_single_enum(k, n, bound)
        return log(N, 2) - 1

    @staticmethod
    @cached_function
    def log_single_enum(n, bound="lower", **kwargs):
        return LinearPruning.log_cost(n, bound="lower", **kwargs)

    # @staticmethod
    # def log_Hk(k, n, bound="lower", pprime=1, **kwargs):
    #     n = ZZ(n)
    #     prob_sol_not_pruned = 1/n
    #     nbases = pprime / prob_sol_not_pruned
    #     return LinearPruning.log_Hk_single_enum(k, n, bound, **kwargs) + log(nbases, 2)

    # @staticmethod
    # def log_cost(n, bound="lower", pprime=1, **kwargs):
    #     """
    #     Sum of expected necessary trees to succeed with probability pprime.
    #     """
    #     n = ZZ(n)
    #     prob_sol_not_pruned = 1/n
    #     nbases = pprime / prob_sol_not_pruned
    #     return LinearPruning.log_single_enum(n, bound) + log(nbases, 2)

    #    log_avg_Skh, log_avg_subtree_N_kD:
    #        1) do basic ones using p 

    @staticmethod
    @cached_function
    def log_avg_Skh(k, h, n, bound="expected", **kwargs):
        """
        Subtree level for single tree, would succeed with probability 1/n.
        """
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        
        if bound == "upper":
            return LinearPruning.log_Hk(k+h, n, "upper") - LinearPruning.log_Hk(k, n, "lower")
        elif bound == "lower":
            return LinearPruning.log_Hk(k+h, n, "lower") - LinearPruning.log_Hk(k, n, "upper")
        elif bound == "expected":
            # for "expected" we'll assume the polynomial speedup happens
            return LinearPruning.log_Hk(k+h, n, "lower") - LinearPruning.log_Hk(k, n, "lower")
        else:
            raise ValueError("`bound` can only be \"lower\", \"upper\" or \"expected\"")

    @staticmethod
    @cached_function
    def log_avg_Skh_single_enum(k, h, n, bound="expected", pprime=1, **kwargs):
        return LinearPruning.log_avg_Skh(k, h, n, bound=bound, **kwargs)

    # @staticmethod
    # def log_avg_Skh(k, h, n, bound="expected", pprime=1, **kwargs):
    #     """
    #     Subtree level for single tree, would succeed with probability 1/n.
    #     """
    #     k, h, n = ZZ(k), ZZ(h), ZZ(n)
    #     prob_sol_not_pruned = 1/n
    #     nbases = pprime / prob_sol_not_pruned
    #     return LinearPruning.log_avg_Skh_single_enum(k, h, n, bound=bound, **kwargs) + log(nbases, 2)

    @staticmethod
    @cached_function
    def avg_children_Ck(k, n, **kwargs):
        return 2**LinearPruning.log_avg_Skh(k, 1, n, **kwargs)

    @staticmethod
    @cached_function
    def log_avg_subtree_N_kD(k, D, n, bound="expected", **kwargs):
        """
        Subtree size for single tree, would succeed with probability 1/n.

        Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
        assumes unit olume
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** LinearPruning.log_avg_Skh(k, h, n, bound=bound)
        return log(N, 2)

    #   log_avg_Skh, log_avg_subtree_N_kD:
    #       2) then implement a "log_cost_rinse_and_repeat" style thing given pprime = m /n and 1)

    # @staticmethod
    # def log_avg_Skh_total(k, h, n, bound="expected", pprime=1, **kwargs):
    #     k, h, n = ZZ(k), ZZ(h), ZZ(n)
    #     prob_sol_not_pruned = 1/n
    #     nbases = pprime / prob_sol_not_pruned
    #     return LinearPruning.log_avg_Skh(k, h, n, bound) + log(nbases, 2)

    # @staticmethod
    # def log_avg_subtree_N_kD_total(k, D, n, bound="expected", pprime=1, **kwargs):
    #     """
    #     Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
    #     assumes unit olume
    #     """
    #     k, D, n = ZZ(k), ZZ(D), ZZ(n)
    #     N = RRRR(0)
    #     for h in range(1, D+1):
    #         N += RRRR(2) ** LinearPruning.log_avg_Skh_total(k, h, n, bound, pprime)
    #     return log(N, 2)


@cached_function
def _alpha_k(alpha, k, n):
    """ Note, this is not GSA's alpha, but rather alpha as in Corollary 2
        and Eq. 16 of [ANSS18]. 
    """
    # from sage.all import RealDistribution
    # alpha, k, n = RRRR(alpha), ZZ(k), ZZ(n)
    # a = k/2
    # b = 1 + (n-k)/2
    # T = RealDistribution('beta', [a, b])

    # #print(f" > > T {T}")

    # T_sample = T.cum_distribution_function_inv(alpha)

    # #print(f" > > T sample {T_sample}")

    # return T_sample
    
    #print(f"alpha {alpha}, k {k}, n {n}")

    return _higher_prec_alpha_k(alpha, k, n)

@cached_function
def _higher_prec_alpha_k(alpha, k, n, max_iterations=64, convergence_err=1e-20):
    """ Note, this is not GSA's alpha, but rather alpha as in Corollary 2
        and Eq. 16 of [ANSS18]. 
    """
    a = k/2
    b = 1 + (n-k)/2

    #alpha = 1/(2**24)

    #print(f" . . (_higher_prec_alpha_k) alpha {alpha}, k {k}, n {n}")

    return _alfk.lib.my_gsl_cdf_beta_Pinv(alpha, a, b, max_iterations, convergence_err)


class CilinderPruningLowerBound:
    """
        References

        [ANSS18] Aono, Y., Nguyen, P.Q., Seito, T., Shikata, J. (2018).
        Lower Bounds on Lattice Enumeration with Extreme Pruning.
        In: Shacham, H., Boldyreva, A. (eds) Advances in Cryptology – CRYPTO 2018.
        Lecture Notes in Computer Science(), vol 10992. Springer, Cham.
        https://doi.org/10.1007/978-3-319-96881-0_21
    """

    @staticmethod
    @cached_function
    def alpha_k(alpha, k, n):
        return _alpha_k(alpha, k, n)

    @staticmethod
    @cached_function
    def log_Hk(k, n, p=None, **kwargs):
        """
        This is Eq (16) of [ANSS18], assuming no constraints on p.

        :param p: lower bound on the single-enum success probability, required unless kwargs contains pprime and nbases, and bool(p) = False
        """
        if k == 0:
            # only one node on root level
            return 0

        if (not p) and ('pprime' in kwargs.keys() and 'nbases' in kwargs.keys()):
            p = kwargs['pprime'] / kwargs['nbases']

        p, k, n = RRRR(p), ZZ(k), ZZ(n)

        #print(f" > > p {p}, k {k}, n {n}")

        alpha_k = CilinderPruningLowerBound.alpha_k(p, k, n)

        #print(f" > > alpha_k {float(alpha_k)}")

        exp_2 = NoPruning.log_Hk(k, n) + k * log(alpha_k, 2) / 2

        #print(f" > > term with alpha_k {float(k * log(alpha_k, 2) / 2)}")
        return float(exp_2)

    @staticmethod
    @cached_function
    def log_Hk_single_enum(k, n, p=None, **kwargs):
        return float(CilinderPruningLowerBound.log_Hk(k, n, p=p, **kwargs))

    @staticmethod
    @cached_function
    def log_cost(n, nbases, pprime, **kwargs):
        # This is the cost for a _single tree_
        n = ZZ(n)
        p = pprime/nbases
        # cost one enumeration
        N = RRRR(0)
        for k in range(1, n+1):
            N += RRRR(2) ** CilinderPruningLowerBound.log_Hk_single_enum(k, n, p, **kwargs)
        # account for repeats
        return float(log(N, 2)-1)

    @staticmethod
    @cached_function
    def log_single_enum(n, nbases, pprime=1, **kwargs):
        return CilinderPruningLowerBound.log_cost(n, nbases, pprime=pprime, **kwargs)

    # @staticmethod
    # def log_Hk(k, n, **kwargs):
    #     k, n = ZZ(k), ZZ(n)
    #     if k == 0:
    #         # only one node on root level
    #         return 0

    #     p = pprime/nbases
    #     return CilinderPruningLowerBound.log_Hk_single_enum(k, n, p) + log(nbases, 2)


    # @staticmethod
    # def log_cost(n, nbases, pprime=1, **kwargs):
    #     """
    #     This computes the overall probability using a "rinse-and-repeat" approach.
    #     The value output here should be lower-bounded by the value output by `log_cost_lb`.

    #     :param pprime: lower bound on the multiple-enum success probability
    #     :param nbases: number of available bases for rinse-and-repeat
    #     """
    #     return CilinderPruningLowerBound.log_single_enum(n, nbases, pprime) + log(nbases, 2)


    # Probably the two _lb that follow are too harsh: they try to do
    # a thing that accounts fo all possible m...
    # May be more meaningful to just choose various values of m instead.

    # The following _lb functions are commented as they try to deal with the multitree automatically,
    # while we have an m-tree above us and wrapper logic, which is therefore incompatible with this
    # def log_Hk_lb(k, n, pprime, **kwargs):
    #     """ This is the "single additive term" of Eq (17)
    #     :param pprime: lower bound on the multiple-enum success probability
    #     """
    #     k, n = ZZ(k), ZZ(n)
    #     exp_2 = NoPruning.log_Hk_single_enum(k, n) + log(pprime/2,2) + log(k,2) + log(beta_function(k/2, 1 + (n-k)/2),2)
    #     return exp_2

    # @staticmethod
    # def log_cost_lb(n, pprime, **kwargs):
    #     """
    #     This is Eq (17) of [ANSS18], assuming we have for free as many bases
    #     as necessary to reach overall probability pprime
    
    #     :param pprime: lower bound on the multiple-enum success probability
    #     """
    #     n = ZZ(n)
    #     N = RRRR(0)
    #     for k in range(1, n+1):
    #         N += RRRR(2) ** CilinderPruningLowerBound.log_Hk_lb(k, n, pprime)
    #     return log(N, 2)-1

    #    log_avg_Skh, log_avg_subtree_N_kD:
    #        1) do basic ones using p. Recommended value: p = (global success prob) / (number of bases available). Eg, p = 1/m for various m = number of bases available.

    # NOTE: we are using a lower bound for log_HK_single_enum since it comes from equation (16) but
    # then we are also using the same style of lower bound for the denominator! Should we use
    # NoPruning instead???

    @staticmethod
    @cached_function
    def log_avg_Skh_single_enum(k, h, n, p=None, bound="upper", **kwargs):
        """
        Subtree level for single tree, would succeed with probability p = the probability of one single enum inside extreme pruning.
        """

        if (not p) and ('pprime' in kwargs.keys() and 'nbases' in kwargs.keys()):
            p = kwargs['pprime'] / kwargs['nbases']

        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        assert(bound == "upper")
        if bound == "lower":
            print(f"CPLB.log_avg_Skh_single_enum(bound=lower, k={k}, h={h}, n={n})")
            return float(CilinderPruningLowerBound.log_Hk_single_enum(k+h, n, p) - NoPruning.log_Hk(k, n))
        elif bound == "upper":
            print(f"CPLB.log_avg_Skh_single_enum(bound=upper, k={k}, h={h}, n={n})")
            return float(CilinderPruningLowerBound.log_Hk_single_enum(k+h, n, p) - CilinderPruningLowerBound.log_Hk_single_enum(k, n, p))
        else:
            raise ValueError(f"`bound` can only be \"lower\" or \"upper\", not \"{bound}\"")

    @staticmethod
    @cached_function
    def log_avg_subtree_N_kD_single_enum(k, D, n, p, **kwargs):
        """
        Subtree size for single tree, would succeed with probability p.

        Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
        assumes unit volume
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        logN = float(-infinity)
        two = RRRR(2)
        for h in range(1, D+1):
            logN = float(log(two**logN + two**CilinderPruningLowerBound.log_avg_Skh_single_enum(k, h, n, p, **kwargs), 2))
        return logN
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** CilinderPruningLowerBound.log_avg_Skh_single_enum(k, h, n, p, **kwargs)
        return float(log(N, 2))

    #   log_avg_Skh, log_avg_subtree_N_kD:
    #       2) then implement a "log_cost_rinse_and_repeat" style thing given pprime = m * p and 1)

    @staticmethod
    @cached_function
    def log_avg_Skh(k, h, n, nbases, pprime, **kwargs):
        """ Skh given n bases and a total success probability pprime. This will compute the single-enum success prob p from pprime and nbases """
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        p = pprime/nbases
        return CilinderPruningLowerBound.log_avg_Skh_single_enum(k, h, n, p, **kwargs) # + log(nbases, 2)

    @staticmethod
    @cached_function
    def avg_children_Ck(k, n, **kwargs):
        return 2**CilinderPruningLowerBound.log_avg_Skh(k, 1, n, **kwargs)

    @staticmethod
    @cached_function
    def log_avg_subtree_N_kD(k, D, n, nbases, pprime, **kwargs):
        """
        Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
        assumes unit volume.

        Starts from overall probability pprime and number of bases nbases, computes the single-enum prob p, and uses that.
        """
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** CilinderPruningLowerBound.log_avg_Skh(k, h, n, nbases, pprime, **kwargs)
        return log(N, 2)

    #   log_avg_Skh, log_avg_subtree_N_kD:
    #       3) then do a lower-bound one using Lemma 6 from the paper

    # Here it combines (16) with the lowerbound for m from (17).
    # This may be both too harsh (trying various values of m may be better)
    # plus it does not really lowerbound unless you choose an upper bound for the denominator
    # and the only one we have is NoPruning.log_Hk_single_enum which is super strict (and on by default)

    # @staticmethod
    # def log_avg_Skh_lb(k, h, n, pprime, use_noprunig_denom=True, **kwargs):
    #     k, h, n = ZZ(k), ZZ(h), ZZ(n)
    #     if use_noprunig_denom:
    #         return CilinderPruningLowerBound.log_Hk_lb(k+h, n, pprime) - NoPruning.log_Hk(k, n)
    #     else:
    #         # print("Warning, this ratio is intended to be a lower bound, but the denomiator being used _is from a lower bound too_")
    #         # print("This means we probably are note getting a lower bound. May want to consider replacing the denominator with NoPRuning.log_Hk(k,n) instead")
    #         return CilinderPruningLowerBound.log_Hk_lb(k+h, n, pprime) - CilinderPruningLowerBound.log_Hk_lb(k, n, pprime)

    # @staticmethod
    # def log_avg_subtree_N_kD_lb(k, D, n, pprime=1, **kwargs):
    #     """
    #     Returns the size of the average subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D
    #     assumes unit olume
    #     """
    #     k, D, n = ZZ(k), ZZ(D), ZZ(n)
    #     N = RRRR(0)
    #     for h in range(1, D+1):
    #         N += RRRR(2) ** CilinderPruningLowerBound.log_avg_Skh_lb(k, h, n, pprime)
    #     return log(N, 2)

    # The following two "non-average" methods need finishing the analysis started on integrating these volumes.

    # Not Implemented
    @staticmethod
    def log_Skh(norm_g, k, h, n):
        """ Assumes unit volume.
        returns the number of nodes on the (k+h)'s enumeration tree level that descend directly from g
        """
        raise NotImplementedError("Not actually implemented. Corollary 2 from the lower-bounds paper may suffice to lower bound the pruning when R_n^2-||g||^2 is used since ||g|| is constant.")
        k, h, n = ZZ(k), ZZ(h), ZZ(n)
        norm_g = RR(norm_g)
        alpha = gsa_alpha(n)
        # assume unit volume lattice
        # R = gaussian_heuristic(RRR(b1norm)**n * RRR(alpha) **((n-1)*n/2), n)
        R = gaussian_heuristic(1, n)
        b1norm = alpha**(-(n-1)/2)
        vol = vol_ball(h, sqrt(R**2 - norm_g**2))
        prod = 1
        for j in range(1, h+1):
            i = n-k-j+1
            bistar = alpha**(i-1) * b1norm
            prod *= bistar
        exp_2 = log(vol,2) - log(prod,2)
        return exp_2

    # Not Implemented
    @staticmethod
    def log_subtree_N_kD(norm_g, k, D, n):
        """
        Returns the size of the subtree spawend from level k, adding nodes on levels k+1, k+2, ..., k+D by a vector on level k of norm norm_g
        assumes unit olume
        """
        raise NotImplementedError()
        k, D, n = ZZ(k), ZZ(D), ZZ(n)
        norm_g = RR(norm_g)
        N = RRRR(0)
        for h in range(1, D+1):
            N += RRRR(2) ** CilinderPruningLowerBound.log_Skh(norm_g, k, h, n)
        return log(N, 2)

    # The following two methods were originally designed to help upper/lower bounding

    # Not Implemented
    @staticmethod
    def log_vol_cilinder_intersection_lb(k, n, Rn, alpha):
        """ Produces a lower bound to the volume of C_{R_1,...,R_k} as
            done in Eq (16)  (using Corollary 2) of [ANSS18].

            :param alpha:   Probability the shortest vector is not pruned by this cilinder intersection.
            :param Rn:      Maximum radius for the cilinder intersection.
        """
        raise NotImplementedError()
        alpha, k, n = RRRR(alpha), ZZ(k), ZZ(n)
        alpha_k = CilinderPruningLowerBound.alpha_k(alpha, k, n)
        return log(vol_ball(k, Rn * sqrt(alpha_k)), 2)

    # Not Implemented
    @staticmethod
    def log_vol_cilinder_intersection_ub(k, Rn):
        """ Produces an upper bound to the volume of C_{R_1,...,R_k} as
            done in Eq (16)  (using Corollary 2) of [ANSS18].

            TODO: Use the Dirichlet volume from [GNR10] rather than the ball.

            :param alpha:   Probability the shortest vector is not pruned by this cilinder intersection.
            :param Rn:      Maximum radius for the cilinder intersection.
        """
        raise NotImplementedError()
        return log(vol_ball(k, Rn) , 2)

try:
    datafile = {
        # 406: "d406_s0.json",
        # 406: "d406_s1.json",
        406: "d406_s2.json", # best radii
        # 406: "d406_s3.json",
        # 622: "d622_s0.json",
        # 622: "d622_s1.json",
        622: "d622_s2.json", # similar values for all, choosing 2 for consistency
        # 622: "d622_s3.json",
        # 623: "d623_s0.json",
        # 623: "d623_s1.json",
        623: "d623_s2.json", # similar values for all, choosing 2 for consistency
        # 623: "d623_s3.json",
        # 872: "d872_s0.json",
        872: "d872_s2.json", # similar values for all, choosing 2 for consistency
        # 873: "d873_s0.json",
        # 873: "d873_s1.json",
        873: "d873_s2.json", # similar values for all, choosing 2 for consistency
        # 873: "d873_s3.json",
    }
    precomputed_upper_bound_data = {}
    for n in datafile:
        with open(f"kyber_pruning_radii/{datafile[n]}") as f:
            precomputed_upper_bound_data[n] = json.load(f)
    class CilinderPruningUpperBound:
        @staticmethod
        @cached_function
        def log_Hk(k, n, p=None, **kwargs):
            if k == 0:
                # only one node on root level
                return 0

            if (not p) and ('pprime' in kwargs.keys() and 'nbases' in kwargs.keys()):
                p = kwargs['pprime'] / kwargs['nbases']

            p, k, n = RRRR(p), ZZ(k), ZZ(n)

            # here you check if we have a file for this dimension and success probability
            if int(log(p,2)) != -64:
                raise ValueError("Upper bounds not computed for probabilities other than 2^-64")

            admitted_dims = [406, 622, 623, 872, 873]
            if n not in admitted_dims:
                raise ValueError(f"Upper bounds computed only for n in {admitted_dims}")

            # print("loading from precomputed table", int(log(p,2)))
            pruning_data = precomputed_upper_bound_data[n]

            if pruning_data["failed"]:
                raise ValueError(f"Data file for n = {n} contains invalid data.")

            # compute using pruning radii
            # Ri = [RRR(r) for r in pruning_data['Ri']]
            # up_logHk = CylinderPruningUB_log_Hk(k, n, Ri, verbose=False)

            # read from fplll (slightly different estimation, essentially identical)
            # fplll_logHk = [RRR(lhk) + 1 for lhk in pruning_data['log_Hk']] # logHk from fplll's estimation
            up_logHk = float(pruning_data['log_Hk'][k-1]) + 1
            return float(up_logHk)
            raise ValueError(f"About to return: {up_logHk}")

            # up_logHk = CylinderPruningUB_log_Hk(k, n, Ri: List[float], verbose=False)
            # return float(up_logHk)


            alpha_k = CilinderPruningLowerBound.alpha_k(p, k, n)
            exp_2 = NoPruning.log_Hk(k, n) + k * log(alpha_k, 2) / 2
            return float(exp_2)
except:
    pass


def tests(PruningClass=NoPruning, **kwargs):
    """
    Note, these are checks for internal consistency of the functions.
    The CilinderPruning numbers do not include the multiplier by the total number of bases in the cost.
    """
    print("===============")
    print("Testing")
    no_pruning = fit_curve(lambda n: PruningClass.log_cost(n, **kwargs))[0]
    print("asymptotic total cost", no_pruning)
    print()

    # no_pruning = fit_curve(lambda n: PruningClass.log_single_enum(n, **kwargs))[0]
    # print("asymptotic single-enum cost", no_pruning)
    # print()

    print("Sanity checks")

    print("trying to fit using N = avg N_{0,n}(0) / 2")
    no_pruning_N0n = fit_curve(lambda n: PruningClass.log_avg_subtree_N_kD(0, n, n, **kwargs)-1)[0]
    print("asymptotic cost", no_pruning_N0n)
    print()

    if PruningClass == NoPruning:
        print("trying to fit using N = N_{0,n}(0) / 2")
        no_pruning_N0n = fit_curve(lambda n: PruningClass.log_subtree_N_kD(0, 0, n, n, **kwargs)-1)[0]
        print("asymptotic cost", no_pruning_N0n)
        print()

        def log_cost_using_N_0n_on_second_level(n):
            n = ZZ(n)
            N = RRRR(0)
            R = gaussian_heuristic(1, n)
            alpha = gsa_alpha(n)
            b1norm = alpha**(-(n-1)/2)
            bnnorm = alpha**(n-1) * b1norm
            floor_ratio = floor(R/bnnorm)
            for coeff in range(-floor_ratio, floor_ratio+1):
                gn_norm = coeff * bnnorm
                N += RRRR(2) ** PruningClass.log_subtree_N_kD(gn_norm, 1, n-1, n, **kwargs)
            return log(N, 2)-1

        print("trying to fit using N = sum_gn N_{1,n-1}(gn) / 2")
        no_pruning_N0n = fit_curve(log_cost_using_N_0n_on_second_level)[0]
        print("asymptotic cost", no_pruning_N0n)
        print()


    def log_cost_using_avg_N_0n_on_second_level(n):
        n = ZZ(n)
        N = RRRR(2) ** PruningClass.log_Hk(1, n, **kwargs) * RRRR(2) ** PruningClass.log_avg_subtree_N_kD(1, n-1, n, **kwargs)
        return log(N, 2)-1

    print("trying to fit using N = H_1 * avg N_{1,n-1} / 2")
    no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_second_level)[0]
    print("asymptotic cost", no_pruning_N0n)
    print()


    def log_cost_using_avg_N_0n_on_second_level(n):
        n = ZZ(n)
        N = RRRR(2) ** PruningClass.log_avg_Skh(0, 1, n, **kwargs) * RRRR(2) ** PruningClass.log_avg_subtree_N_kD(1, n-1, n, **kwargs)
        return log(N, 2)-1

    print("trying to fit using N = avg S_{0,1} * avg N_{1,n-1} / 2")
    no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_second_level)[0]
    print("asymptotic cost", no_pruning_N0n)
    print()



    lvl = 20
    def log_cost_using_avg_N_0n_on_given_level(n):
        n = ZZ(n)
        N = RRRR(2) ** PruningClass.log_avg_subtree_N_kD(0, lvl, n, **kwargs) + RRRR(2) ** PruningClass.log_avg_Skh(0, lvl, n, **kwargs) * RRRR(2) ** PruningClass.log_avg_subtree_N_kD(0 + lvl, n-lvl, n, **kwargs)
        return log(N, 2)-1

    print("trying to fit using N = avg N_{0, lvl} + avg S_{0,lvl} * avg N_{lvl,n-lvl} / 2")
    no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_given_level)[0]
    print("asymptotic cost", no_pruning_N0n)
    print()


    print("some subtree tests, these may differ in asymptotic to the total cost, but still coincide")

    k = 20
    h = 50
    D = 80
    print("avg_N_{k,D}")
    no_pruning_N0n = fit_curve(lambda n: PruningClass.log_avg_subtree_N_kD(k, D, n, **kwargs))[0]
    print("asymptotic cost", no_pruning_N0n)
    print()

    def log_cost_using_avg_N_0n_on_given_level(n):
        n = ZZ(n)
        N = RRRR(2) ** PruningClass.log_avg_subtree_N_kD(k, h, n, **kwargs) + RRRR(2) ** PruningClass.log_avg_Skh(k, h, n, **kwargs) * RRRR(2) ** PruningClass.log_avg_subtree_N_kD(k + h, D-h, n, **kwargs)
        return log(N, 2)

    print("trying to fit using avg_N_{k,D} = avg N_{k, h} + avg S_{k,h} * avg N_{h,D-h}")
    no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_given_level)[0]
    print("asymptotic cost", no_pruning_N0n)
    print()


    from sage.all import beta as beta_function
    print("Checking children can be predicted")
    print("Number-of-children tests are currently not passing, likely cause of the analysis not using beta() for pruned versions")
    n = 30
    R = gaussian_heuristic(1, n)
    alpha = gsa_alpha(n)
    b1norm = alpha**(-(n-1)/2)
    for k in range(15,20):
        b_n_mk_norm = alpha**(n-k-1) * b1norm
        print(f"On level {k}, (this one may go wrong for PruningClass != NoPruning)")
        # print("avg_Skh_lb_ub(k,1,n) = %.4f" % (2**NoPruning.log_avg_Skh_ub(k, 1, n)))
        # print("W_{k+1}(avg g)}      = %.4f" % (2 * R * sqrt(2/(k+2)) / b_n_mk_norm))
        print("avg W_{k+1}(g)}      = %.4f" % (R * beta_function(k/2 + 1, 1/2) / b_n_mk_norm))
        print("avg Skh_lb(k,1,n)    = %.4f" % (2**PruningClass.log_avg_Skh(k, 1, n, **kwargs)))
        print("avg Ck(k,1,n)        = %.4f" % (PruningClass.avg_children_Ck(k, n, **kwargs)))
        print()


def linear_pruning_plots(ymax=None):
    from sage.all import line, save
    print("===============")
    print("Plotting Linear Pruning")
    print()

    f = lambda x, exponent: x * log(x,2) * exponent['nl2'] + x * exponent['n'] + log(x,2) * exponent['l2'] + exponent['c']
    g = line([], title=f"Linear pruning", ymin=0, ymax=ymax)

    no_pruning = fit_curve(NoPruning.log_cost)[0]
    g += line([(x, f(x, exponent=no_pruning)) for x in range(50, 1000)], linestyle="--", legend_label="no pruning", color='red', title=f"Linear pruning", ymin=0, ymax=ymax)
    print(f"asymptotic cost no-pruning", no_pruning)
    print()

    for bound in ["upper", "lower"]:
        colors = (x for x in ['orange', 'green', 'purple', 'black', 'brown', 'violet'])

        print(f"Linear pruning (heuristic {bound} bound)")
        lin_pruning_single = fit_curve(lambda n: LinearPruning.log_single_enum(n, bound))[0]
        g += line([(x, f(x, exponent=lin_pruning_single)) for x in range(50, 1000)], linestyle="solid" if bound == "upper" else "--", legend_label=f"single enum ({bound}, bound)", color='blue', ymin=0, ymax=ymax, title=f"Linear pruning")
        print("Single enum", lin_pruning_single)
        print("Total cost")
        for pprime in [x / 100 for x in range(1, 100, 30)]:
            pprime = .61 ################################
            lin_pruning = fit_curve(lambda n: LinearPruning.log_cost(n, bound, pprime))[0]
            g += line([(x, f(x, exponent=lin_pruning)) for x in range(50, 1000)], linestyle="solid" if bound == "upper" else "--", legend_label=f"p = {pprime} enum ({bound}, bound)", color=next(colors), ymin=0, ymax=ymax, title=f"Linear pruning")
            print(f"  p_succ = {pprime},", lin_pruning)
            print()
            break ################################
        print()

    save(g, "linear_prun.png", dpi=200)


def cilinder_pruning_plots(ymax=None):
    from sage.all import line, save
    print("===============")
    print("Plotting Cilinder Pruning")

    print("using lower bound to rinse-and-repeat")
    for pprime in [x / 100 for x in range(1, 100, 30)]:
        # pprime = 0.61 ####################################
        colors = (x for x in ['red', 'blue', 'orange', 'green', 'purple', 'black',
                          'brown', 'violet'])
        g = line([], title=f"extreme cilinder pruning lower bound, p' = {pprime}", ymin=0, ymax=ymax)
        f = lambda x, exponent: x * log(x,2) * exponent['nl2'] + x * exponent['n'] + log(x,2) * exponent['l2'] + exponent['c']

        no_pruning = fit_curve(NoPruning.log_cost)[0]
        print(f"asymptotic cost no-pruning", no_pruning)
        print()
        g += line([(x, f(x, exponent=no_pruning)) for x in range(50, 1000)], linestyle="--", legend_label="no pruning", color='red', ymin=0, ymax=ymax, title=f"extreme cilinder pruning lower bound, p' = {pprime}")



        # log_cost_lb = lambda x: CilinderPruningLowerBound.log_cost_lb(x, pprime=pprime)
        # no_pruning = fit_curve(log_cost_lb)[0]
        # print(f"asymptotic cost lb (p' = {pprime})", no_pruning)
        # print()
        # g += line([(x, f(x, exponent=no_pruning)) for x in range(50, 1000)], linestyle="--", legend_label="lb", ymin=0, ymax=ymax, title=f"extreme cilinder pruning lower bound, p' = {pprime}")


        print("assuming explicit rinse and repeat")
        for log10_nbases in [0, 1, 2, 3, 5, 10, 20]:
            # log10_nbases = 10 ####################################
            nbases = 10**log10_nbases
            log_cost_m = lambda x: CilinderPruningLowerBound.log_cost(x, nbases, pprime=pprime) + log(nbases, 2)
            no_pruning = fit_curve(log_cost_m)[0]
            print(f"asymptotic cost, (log10(nbases), p') = ({log10_nbases}, {pprime})", no_pruning)
            _color = next(colors)
            g += line([(x, f(x, exponent=no_pruning)) for x in range(50, 1000)], legend_label=f"log10(m) = {log10_nbases}", color=_color, ymin=0, ymax=ymax, title=f"extreme cilinder pruning lower bound, p' = {pprime}")

            single_red_cost = lambda x: CilinderPruningLowerBound.log_single_enum(x, nbases, pprime=pprime)
            no_pruning = fit_curve(single_red_cost)[0]
            print(f"asymptotic cost, (log10(nbases), p'), one basis = ({log10_nbases}, {pprime})", no_pruning)
            print()
            g += line([(x, f(x, exponent=no_pruning)) for x in range(50, 1000)], linestyle="dotted", legend_label=f"log10(m) = {log10_nbases}, 1 base", color=_color, ymin=0, ymax=ymax, title=f"extreme cilinder pruning lower bound, p' = {pprime}")

            # print("computing cost using average subtree size from second level")
            # print(f"pprime = {pprime}, nbases = {nbases}")
            # print()
            # print("doing rinse and repeat")
            # def log_cost_using_avg_N_0n_on_second_level(n):
            #     n = ZZ(n)
            #     N = RRRR(2) ** CilinderPruningLowerBound.log_Hk(1, n, pprime/nbases) * RRRR(2) ** CilinderPruningLowerBound.log_avg_subtree_N_kD_total(1, n-1, n, nbases, pprime)
            #     return log(N, 2)-1

            # no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_second_level)[0]
            # print("asymptotic cost", no_pruning_N0n)
            # print()
            # print("using (17) lower bound, no-pruning denominator")
            # def log_cost_using_avg_N_0n_on_second_level(n):
            #     n = ZZ(n)
            #     N = RRRR(2) ** CilinderPruningLowerBound.log_Hk_lb(1, n, pprime) * RRRR(2) ** CilinderPruningLowerBound.log_avg_Skh_lb(1, n-1, n, pprime, use_noprunig_denom=True)
            #     return log(N, 2)-1

            # no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_second_level)[0]
            # print("asymptotic cost", no_pruning_N0n)
            # print()
            # print("using (17) lower bound, lb denominator")
            # def log_cost_using_avg_N_0n_on_second_level(n):
            #     n = ZZ(n)
            #     N = RRRR(2) ** CilinderPruningLowerBound.log_Hk_lb(1, n, pprime) * RRRR(2) ** CilinderPruningLowerBound.log_avg_Skh_lb(1, n-1, n, pprime, use_noprunig_denom=False)
            #     return log(N, 2)-1

            # no_pruning_N0n = fit_curve(log_cost_using_avg_N_0n_on_second_level)[0]
            # print("asymptotic cost", no_pruning_N0n)
            # print()
            # break ####################################

        print()
        print()

        save(g, f"pprime_{pprime}.png", dpi=200)
        # break ####################################


def plot_children():
    from sage.all import line, save
    colors = (x for x in ['orange', 'green', 'purple', 'black', 'brown', 'violet'] * 3)
    n = 400
    g = line([], title=f"Log2 Avg Children(k), n = {n}")
    for PruningClass in [NoPruning, LinearPruning]:
        g += line(
            [(k, log(PruningClass.avg_children_Ck(k, n),2)) for k in range(n)],
            color=next(colors),
            legend_label=f"{PruningClass.__name__}",
        )
    for pprime in [.1, .5, .9]:
        for nbases in [1, 100, 10**10, 10**20]:
            g += line(
                [(k, log(CilinderPruningLowerBound.avg_children_Ck(k, n, pprime=pprime, nbases=nbases),2)) for k in range(n)],
                color=next(colors),
                legend_label=f"CilinderPLB, p'= {pprime}, lg m = %.2f" % log(nbases,2),
            )
    save(g, "avg_children.png", dpi=200)


# missing="""
# A class that helps dealing with the top m-tree
# """
# print(missing)


def test_gap_in_lb_up_denom(PruningClass, n=100, k=20, h=50, D=80, **kwargs):
    bound = kwargs['bound']
    print("log_avg_S_{k,h} (%s) = %.2f" % (bound, PruningClass.log_avg_Skh(k, h, n, **kwargs)))
    print("log_avg_N_{k,D} (%s) = %.2f" % (bound, PruningClass.log_avg_subtree_N_kD(k, h, n, **kwargs)))



def bound_plots(test):
    try:
        dim, prec, nbases = test
        # print(get_pruning_bound(20, 200, .41))

        """
        NOTE: log_Hk assumes the GH from GSA ratio. We are passing Ri as absolute values. We have to turn them into scales, and integrate them with log_Hk computations.
        """

        target_pr = float(1/nbases)
        if nbases == 1:
            target_pr = 0.999

        fn = f"{dim}-{prec}-{nbases}.json"
        if pathlib.Path(fn).exists():
            with open(fn) as f:
                pruning_data = json.load(f)
        else:
            pruning_data = get_pruning_bound(dim, prec, target_pr)
            with open(fn, 'w') as f:
                json.dump(pruning_data, f)
    except:
        with open(f"kyber_pruning_radii/{test[0]}") as f:
            pruning_data = json.load(f)
        dim = pruning_data["dimension"]
        prec = pruning_data["precision"]
        nbases = round(1/float(pruning_data["probability"]))
        target_pr = float(pruning_data["probability"])
        fn = f"{test[0]}.json"

    if pruning_data['failed']:
        print("Something went wrong. Low precision?")
        exit(1)
    if "log_Hk" not in pruning_data:
        return

    R = gaussian_heuristic(1, dim)
    Ri = [RRR(r) for r in pruning_data['Ri']]
    logHk = [RRR(lhk) + 1 - 0*float(log(nbases,2)) for lhk in pruning_data['log_Hk']]

    fplll_single_cost = RRR(pruning_data["log_cost"]) + 1 - 0*float(log(nbases,2))
    fplll_prob = RRR(pruning_data["probability"])

#   print(float(2**fplll_single_cost))
#   print(float(sum(2**lhk for lhk in logHk)))
#   print(abs(2**fplll_single_cost - sum(2**lhk for lhk in logHk)))
    assert(abs(2**fplll_single_cost - sum(2**lhk for lhk in logHk))/float(2**fplll_single_cost) < 0.01)

    # V = 1
    # R = gaussian_heuristic(1, dim)
    # # Ri = [R * .1, R * .2, R * .3, R * .4, R * .5, R * .6, R * .7, R * .8, R * .9, R]
    # Ri = list(float(R * (_/dim)**1) for _ in range(1, dim+1)) # linear pruning

    print()
    print("NoPruning.log_cost(dim, R=Ri[-1])")
    print()
    print(NoPruning.log_cost(dim, R=Ri[-1]).n())

    print()
    print("CilinderPruningUpperBound.log_cost(dim, str(Ri))")
    print()
    print("fplll single enum:", float(fplll_single_cost))
    print("fplll prob :", float(fplll_prob))
    print("fplll total:", fplll_single_cost - log(fplll_prob, 2).n())
    print()

    our_prob = Pr_Chen13_succ_prob_exact(Ri)
    our_single, our_total = CylinderPruningUB_log_cost(dim, Ri)
    our_single_sum = log(sum(2**CylinderPruningUB_log_Hk(k, dim, Ri) for k in range(1, dim+1)), 2).n() - 1
#   our_single_sum = log(sum(2**CylinderPruningUB_log_Hk(k, dim, [gaussian_heuristic(1, dim) * Ri[i] / Ri[dim-1] for i in range(dim)], normalise_Ri=False) for k in range(1, dim+1)), 2).n() - 1
    print("our single:", float(our_single))
    print("our single:", float(our_single_sum))
    print("our prob   :", our_prob)
    print("our prob fp:", Pr_Chen13_succ_prob_exact(Ri, fplll=True))
    print("our total:", our_single - log(our_prob, 2).n())
    print("our total:", our_total)
    print()

    print()
    print("CilinderPruningLowerBound.log_cost(dim, pprime=1, nbases=10)")
    print()
    single_enum_lb = CilinderPruningLowerBound.log_cost(dim, pprime=1/nbases, nbases=1)
    total_enum_lb = CilinderPruningLowerBound.log_cost(dim, pprime=1, nbases=nbases)
    print("single enum:", single_enum_lb)
    print("repeat:", float(log(nbases, 2)))
    print("total: ", total_enum_lb)


    from sage.all import line, save

    fp_pts = []
    our_pts = []
    lb_pts = []
    np_pts = []

    for k in range(1, dim+1, 1):
        fplll_lhk = float(logHk[k-1])
        our_lhk = float(CylinderPruningUB_log_Hk(k, dim, Ri, verbose=False)) #- float(log(Pr_Chen13_succ_prob_exact(Ri), 2))
        lb_lhk = CilinderPruningLowerBound.log_Hk(k, dim, pprime=1/nbases, nbases=1)
        np_lhk = NoPruning.log_Hk(k, dim)
        print("k", k, "\t", fplll_lhk, "\t", our_lhk, "\t", lb_lhk)
        if not math.isinf(fplll_lhk):
            fp_pts += [(k, fplll_lhk)]
        our_pts += [(k, our_lhk)]
        lb_pts += [(k, lb_lhk)]
        np_pts += [(k, np_lhk)]

    print(float(log(sum(2**pts[1] for pts in lb_pts),2))-1)
    print(single_enum_lb)
    assert(abs(log(sum(2**pts[1] for pts in lb_pts),2)-1 - single_enum_lb) < 0.0001)

    g = line([])
    g += line(fp_pts, color="blue", legend_label="fplll")
    g += line(our_pts, color="green", legend_label="our")
    g += line(lb_pts, color="red", legend_label="lb")
    g += line(np_pts, color="black", legend_label="np")
    save(g, f"{fn}.png")

    gg = line([], title="square radii")
    gg += line([(k, Ri[k-1]**2) for k in range(1,dim+1)], color="blue", legend_label="fplll")
    gg += line([(k, CilinderPruningLowerBound.alpha_k(target_pr, k, dim)) for k in range(1,dim+1)], color="red", legend_label="lb")
    gg += line([(1, 1), (dim, 1)], color="black", legend_label="np")
    save(gg, f"radii-{fn}.png")

    return


import json
def cylinder_pruning_ub_tests(eps=2**-20, biteps=4):

    tests = [
    # (48, 220, 1),
    # (48, 220, 1/0.6),
    # (48, 220, 1/0.999),
    # (48, 220, 2),
    # (48, 220, 4),
    # (48, 220, 2**16),
    # (120, 300, 1),
    # (120, 300, 1/0.6),
    # (120, 300, 1/0.999),
    # (120, 300, 2),
    # (120, 300, 4),
    # (120, 300, 2**16),
    # (120, 300, 2**32),
    # (160, 300, 1),
    # (160, 300, 1/0.6),
    # (160, 300, 1/0.999),
    # (160, 300, 2),
    # (160, 300, 4),
    # (160, 300, 2**16),
    # (160, 300, 2**32),
        ("d406_s0.json",),
        ("d406_s1.json",),
        ("d406_s2.json",),
        ("d406_s3.json",),
        ("d622_s0.json",),
        ("d622_s1.json",),
        ("d622_s2.json",),
        ("d622_s3.json",),
        ("d623_s0.json",),
        ("d623_s1.json",),
        ("d623_s2.json",),
        ("d623_s3.json",),
        ("d872_s0.json",),
        ("d872_s2.json",),
        ("d873_s0.json",),
        ("d873_s1.json",),
        ("d873_s2.json",),
        ("d873_s3.json",),
    ]

    for t in tests:
        bound_plots(t)

    return
    datadir = "kyber_pruning_radii"
    profiles = [
        # "example_high_avg.json", # 1/2 prob, odd k averages Hk
        # "example_high_logs.json", # 1/2 prob, odd k averages log_Hk
        # "example_low_avg.json", # low prob, odd k averages Hk
        # "example_low_logs.json", # low prob, odd k averages log_Hk
        # "150_1010.json",
        # "150_6827.json",
        # "150_1010_2.json",
        # "150_6827_2.json",
        # "160_1010_2.json",
        # "160_6827.json",
        # "160_6827_samepr.json",
        "160_1010_samepr.json",
        # "d406_s0.json",
        # "d406_s1.json",
        # "d406_s2.json",
        # "d406_s3.json",
        # "d623_s0.json",
        # "d623_s1.json",
        # "d623_s2.json",
        # "d873_s0.json"
    ]
    for profile in profiles:
        print("checking", profile)
        with open(f"{datadir}/{profile}") as f:
            pruning_data = json.load(f)
        print(list(pruning_data.keys()))
        if pruning_data["failed"]:
            raise Exception(f"Pruning data for {profile} failed")
        dim = pruning_data["dimension"]
        vol = RRR(pruning_data["volume"])
        Ri = [RRR(_) for _ in pruning_data["Ri"]]
        single_enum_pr = RRR(pruning_data["probability"])
        fplll_single_enum_cost = RRR(pruning_data["log_cost"])
        pruner = pruning_data["pruner"]
        R = gaussian_heuristic(vol, dim)
        nbases = round(float(-log(single_enum_pr, 2)))
        fplll_log_Hk = [RRR(_)+1. for _ in pruning_data["log_Hk"]] # undo /2 due to symmetry
        print("nbases", nbases)
        # assert(nbases == 64)

        # print("dim", dim)
        # print("vol", vol)
        # print("single_enum_pr", single_enum_pr)
        # print("pruner", pruner)
        # print("log2 nbases", log(1/single_enum_pr,2).n())

        # sanity checks

        # fplll prob and our prob match
        our_single_enum_pr = CilinderPruningUpperBound.Pr_Chen13_succ_prob_exact(Ri)
        assert(abs(our_single_enum_pr - single_enum_pr) < eps)
        print("win prob check passed")

        # print("sum", float(log(sum(2**(lhk) for lhk in fplll_log_Hk), 2)))
        # print("sum", float(fplll_single_enum_cost+1))
        assert(abs(log(sum(2**(lhk) for lhk in fplll_log_Hk), 2) - (fplll_single_enum_cost+1)) < eps)

        # print("our tot", CilinderPruningUpperBound.log_cost(dim, Ri))
        # exit(0)

        from sage.all import line, save
        g = line([])
        g += line([(i, fplll_log_Hk[i]) for i in range(dim)], color="blue", legend_label="fplll")
        g += line([(i, CilinderPruningLowerBound.log_Hk(i+1, dim, RRR(1/nbases))) for i in range(dim)], color="red", legend_label="lb")
        g += line([(i, CilinderPruningUpperBound.log_Hk(i+1, dim, Ri, normalise_Ri=False)) for i in range(dim)], color="green", legend_label="ub")
        save(g, f"{profile}.png", dpi=150)
        continue

        # our LB costs <= our UB costs <= our nopruning costs
        falses = 0
        for k in range(2, dim+1, 2):
            ub_cost_k = CilinderPruningUpperBound.log_Hk(k, dim, Ri)
            # lb_cost_k = CilinderPruningLowerBound.log_Hk(k, dim, RRR(1/nbases))
            # print(float(ub_cost_k))
            # print(float(lb_cost_k))
            # print(float(fplll_log_Hk[k-1]))
            print(ub_cost_k)
            # input()
            # print(float(fplll_log_Hk[k-1]) > float(lb_cost_k))
            # if float(fplll_log_Hk[k-1]) < float(lb_cost_k):
            #     falses += 1
            # input()
            # np_cost_k = NoPruning.log_Hk(k, dim)
            # if lb_cost_k > ub_cost_k:
            #     print("comp with lb fails")
            #     print("lb", RRR(lb_cost_k))
            #     print("ub", RRR(ub_cost_k))
            #     input()
            # try:
            #     if ub_cost_k > np_cost_k:
            #         print("comp with np fails")
            # except Exception as e:
            #     print("ub_cost_k", ub_cost_k.n())
            #     print("np_cost_k", np_cost_k.n())
            #     raise e
            # print(f"k = {k} done")
        # print("falses:", falses)
        exit(0)
        ub_cost_total = CilinderPruningUpperBound.log_cost(dim, Ri)
        lb_cost_total = CilinderPruningLowerBound.log_cost(dim, nbases, 1.)
        np_cost_total = NoPruning.log_cost(dim)
        assert(lb_cost_total <= ub_cost_total)
        assert(ub_cost_total <= np_cost_total)

        # fplll cost and our costs are relatively close
        our_single_cost = CilinderPruningUpperBound.log_cost(dim, Ri, single_enum=True)
        print("our_single_cost", our_single_cost)
        print("fplll_single_cost", fplll_single_enum_cost)
        assert(fplll_single_enum_cost >= our_single_cost)
        assert(abs(fplll_single_enum_cost - our_single_cost) < biteps)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-beta', type=int, default=30, help="BKZ block size")
    args = parser.parse_args()

    # print("\n\nNo Pruning\n")
    # tests(NoPruning)
    # # print(f"\n\nLinearPruning\n")
    # # tests(LinearPruning)
    # pprime = .9
    # nbases = 10**10
    # use_noprunig_denom = False
    # print(f"\n\nCilinderPruningLowerBound, p' = {pprime}, m = {nbases}, use_nopruning_denom = {use_noprunig_denom}\n")
    # tests(CilinderPruningLowerBound, pprime=pprime, nbases=nbases, use_noprunig_denom=use_noprunig_denom)

    cylinder_pruning_ub_tests()
    exit(0)

    # linear_pruning_plots(ymax=800)
    cilinder_pruning_plots(ymax=800)
    exit(0)

    for n in [100, 400, 800]:
        D = round(n/4)
        h = round(D * 2 / 3)
        for _k in range(4):
            k = round(_k * n / 4)
            for bound in ["lower", "upper", "expected"]:
                print(f"\nLinearPruning, bound={bound}, (n, k, h, D) = ({n}, {k}, {h}, {D})")
                test_gap_in_lb_up_denom(LinearPruning, n=n, k=k, h=h, D=D, bound=bound)
            print("-------------------")
 
    for n in [100, 400, 800]:
        D = round(n/4)
        h = round(D * 2 / 3)
        for _k in range(4):
            k = round(_k * n / 4)
            for bound in ["lower", "upper"]:
                print(f"\nCilinderPruningLowerBound, bound={bound}, (n, k, h, D) = ({n}, {k}, {h}, {D})")
                pprime = .5
                nbases = 10**4
                test_gap_in_lb_up_denom(CilinderPruningLowerBound, n=n, k=k, h=h, D=D, bound=bound, pprime=pprime, nbases=nbases)
            print("-------------------")
