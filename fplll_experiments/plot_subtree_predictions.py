import os 
import json
import profile
from random import gauss
from matplotlib import legend_handler
from sage.all import line, save, sqrt, log, prod, beta, list_plot, multi_graphics
from bkz_simulators.sim import Sim, LLLProfile
import sys
sys.path.append("../LatticeHeuristics")
import enumTrees

def gaussian_heuristic(v, n, verbose=True):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """
    from math import pi, sqrt, e
    return ((pi * n)**(1./(2.*n))) * sqrt(n/(2 * pi * e)) * (v ** (1./n))


def lwe_gh(lwe):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    vol = q**m
    dim = n+m+1
    return gaussian_heuristic(vol, dim)


def profile_gh(profile, verbose=False):
    # profile contains square norms
    dim = len(profile)
    vol = sqrt(prod(profile))
    if verbose:
        print("dim", dim)
        print("vol", vol)
        print("gh", gaussian_heuristic(vol, dim))
    return gaussian_heuristic(vol, dim)


def test_discrepancy(
        full_sim_profile,
        full_gs2_profile,
        full_dim,
        block_dim,
        block_index, # first index is 0
        logHk,
        CKh,
        gh_factor,
        used_GH
        ):

    def vol_n_ball(n, R):
        from sage.all import gamma, pi, log
        # return n*log(R, 2) +  (n/2)*log(pi,2) - log(gamma(n/2+1),2)
        return float(R**n * pi**(n/2) / gamma(n/2+1))

    def rad_n_ball(n, v):
        from sage.all import pi, gamma, sqrt
        return float(gamma(n/2+1)**(1/n) * v**(1/n) / sqrt(pi))

    print("\n ==== Testing discrepancies:\n")
    print("full dim:", full_dim)
    print("block dim:", block_dim)
    print(f"block indices (basis from 0): [{block_index}, {block_index+block_dim-1}]")
    sim_block_basis = full_sim_profile[block_index:min(block_index+block_dim, full_dim)]
    gs2_block_basis = full_gs2_profile[block_index:min(block_index+block_dim, full_dim)]
    print(f"block basis (sim): [{sim_block_basis[0]}, {sim_block_basis[1]}, ..., {sim_block_basis[-1]}], len = {len(sim_block_basis)}")
    print(f"block basis (gs2): [{gs2_block_basis[0]}, {gs2_block_basis[1]}, ..., {gs2_block_basis[-1]}], len = {len(gs2_block_basis)}")
    from sage.all import log
    l2 = lambda x: float(log(x,2))
    log_vol_sim_block_basis = sum(map(l2, sim_block_basis))/2
    log_vol_gs2_block_basis = sum(map(l2, gs2_block_basis))/2
    print("logvol block basis (sim):", log_vol_sim_block_basis)
    print("logvol block basis (gs2):", log_vol_gs2_block_basis)
    GH_sim_basis = profile_gh(sim_block_basis)
    GH_gs2_basis = profile_gh(gs2_block_basis)
    test_GH_sim_basis = rad_n_ball(block_dim, 2**log_vol_sim_block_basis)
    test_GH_gs2_basis = rad_n_ball(block_dim, 2**log_vol_gs2_block_basis)
    print("GH block (sim):", GH_sim_basis, f"(test: {test_GH_sim_basis})")
    print("GH block (gs2):", GH_gs2_basis, f"(test: {test_GH_gs2_basis})")
    print("used GH:", used_GH)
    print("GH factor:", gh_factor)
    print("sqrt(GH factor):", sqrt(gh_factor))
    print("sim_R", used_GH * sqrt(gh_factor))

    k = block_dim//2
    print("k:", k)
    log_vol_proj_sublattice_sim = sum(map(l2, sim_block_basis[block_dim-k:]))/2
    log_vol_proj_sublattice_gs2 = sum(map(l2, gs2_block_basis[block_dim-k:]))/2
    print("logvol proj sublattice basis (sim):", log_vol_proj_sublattice_sim)
    print("logvol proj sublattice basis (gs2):", log_vol_proj_sublattice_gs2)
    log_Hk_sim = l2(vol_n_ball(k, GH_sim_basis)) - log_vol_proj_sublattice_sim
    log_Hk_gs2 = l2(vol_n_ball(k, GH_gs2_basis)) - log_vol_proj_sublattice_gs2
    print("log_Hk_sim:", log_Hk_sim)
    print("log_Hk_gs2:", log_Hk_gs2)

    print("(fun) logHk (sim)", float(logHk(k, "sim", sim_R=GH_sim_basis)))
    print("(fun) logHk (gs2)", float(logHk(k, "gs2", sim_R=GH_gs2_basis)))
    print("(fun) logHk (gsa)", float(logHk(k, "gsa", sim_R=None)))
    print("(fun) logHk(sim_R) (gsa)", float(logHk(k, "gsa", sim_R=GH_sim_basis)))
    print("\n ==== end of tests\n")


def gen_plots(data, plots_dir, lwe, bkz, subtree_root_level, _tour=None, _block=None, treedata=None, dpi=150):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    sim = Sim(LLLProfile(n, q, m))

    if _tour:
        data_tuple = [data[int(_tour)]]
        tour = int(_tour)
        sim(bkz['beta'], tour)
    else:
        data_tuple = data
        tour = 0

    for tour_stats in data_tuple:
        if _block:
            tour_stats_tuple = [_block]
        else:
            tour_stats_tuple = tour_stats

        if int(tour) in [9, 14]:
            for index in tour_stats_tuple:
                # if int(index) != 100:
                #     continue
                print(f"tour: {tour}, index: {index}, beta: {bkz['beta']}")

                # sim_block_prof = sim.profile[int(index):min(int(index)+bkz['beta'], len(sim.profile))]
                sim_block_prof = sim.profile[int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
                sim_gh = profile_gh(sim_block_prof)
                if treedata:
                    # gs2 = treedata[tour][index]["gs"][int(index):min(int(index)+bkz['beta'], len(sim.profile))]
                    gs2 = treedata[tour][index]["gs"][int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
                sim_R = sqrt(bkz['ghfactor']) * sim_gh
                stats = tour_stats[index]
                avg = stats["avg"][subtree_root_level:]
                var = stats["var"][subtree_root_level:]

                from sage.all import pi, gamma
                k = subtree_root_level
                n = bkz['beta']
                def logHk(k, kind, sim_R=sim_R):
                    if kind == "gsa":
                        if sim_R:
                            sim_V = prod(sim_block_prof)
                            sim_V = sqrt(sim_V)
                        else:
                            sim_V = None
                        return enumTrees.NoPruning.log_Hk(k, bkz['beta'], R=sim_R, V=sim_V)
                    elif kind == "sim":
                        res = log(sim_R**k * pi**(k/2) / gamma(k/2+1),2)
                        covol = prod([sim_block_prof[i-1] for i in range(bkz['beta']-k+1, bkz['beta']+1)])
                        covol = sqrt(covol)
                        res = res - log(covol,2)
                        return res
                    else:
                        res = log(sim_R**k * pi**(k/2) / gamma(k/2+1),2)
                        covol = prod([gs2[i-1] for i in range(bkz['beta']-k+1, bkz['beta']+1)])
                        covol = sqrt(covol)
                        res = res - log(covol,2)
                        return res

                def Ckh(k, h, sim):
                    # formula obtained by direct integration.
                    res = pi**(h/2) * sim_R**h * gamma(k/2+1)/gamma((k+h)/2+1)
                    if sim:
                        covol = prod([sim_block_prof[i-1] for i in range(bkz['beta']-k, bkz['beta']-k-h, -1)])
                    else:
                        covol = prod([gs2[i-1] for i in range(bkz['beta']-k, bkz['beta']-k-h, -1)])
                    covol = sqrt(covol)
                    res = res / covol
                    return res

                try:
                    print(index, k)
                    g = line([])
                    g += line([(k + i + 1, 2**(logHk(k+i+1, "gsa")-logHk(k, "gsa"))) for i in range(len(avg))], xmin=subtree_root_level, linestyle="dashed", color="black", legend_label="$\\mathbb{E}[H_{k+h}]/\\mathbb{E}[H_k]$ gsa")
                    g += line([(k + i + 1, 2**(logHk(k+i+1, "sim")-logHk(k, "sim"))) for i in range(len(avg))], xmin=subtree_root_level, linestyle="dashed", color="green", legend_label="$\\mathbb{E}[H_{k+h}]/\\mathbb{E}[H_k]$ sim")
                    g += line([(k + i + 1, 2**(logHk(k+i+1, False)-logHk(k, False))) for i in range(len(avg))], xmin=subtree_root_level, linestyle="dashed", color="blue", legend_label="$\\mathbb{E}[H_{k+h}]/\\mathbb{E}[H_k]$ gs")

                    g += list_plot([(k + i + 1, avg[i]) for i in range(len(avg))], xmin=subtree_root_level,  marker="D", size=20, color="blue", legend_label="measured avg. $S_{k,h}$", axes_labels=["$k+h$", ""])

                    g.set_legend_options(loc='upper right')

                    if treedata:
                        gs_line = line([(i, log(gs2[i],2)) for i in range(len(gs2))], xmin=1, transparent=False, frame=True, ticks=[[], []], color="blue") # , legend_label="$\log_2(\|b^*_i\|^2$)"
                        gs_line += line([(i, log(sim_block_prof[i],2)) for i in range(len(gs2))], xmin=1, linestyle="dashed", color="green", transparent=False, frame=True, ticks=[[], []]) # , legend_label="[CN11]"
                        g = multi_graphics([g, (gs_line, (0.725, 0.225, 0.2, 0.2))])

                    filename = f"{plots_dir}/{tour}-{index}"
                    print(filename)
                    # save(g, filename+".png", dpi=dpi)
                    save(g, filename+".pdf")
                except:
                    pass
                    # may fail at the last few blocks

        sim(bkz['beta'], 1)
        tour += 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-tree-stats-filename', type=str, default=None)
    parser.add_argument('-subtree-stats-filename', type=str, default=None)
    parser.add_argument('-plot-dir', type=str, default="./plots")
    parser.add_argument('-n', type=int, default=72, help="LWE secret dimension")
    parser.add_argument('-q', type=int, default=97, help="LWE modulo")
    parser.add_argument('-sd', type=float, default=1., help="LWE discrete gaussian standard deviation")
    parser.add_argument('-m', type=int, default=87, help="LWE samples")
    parser.add_argument('-beta', type=int, default=30, help="BKZ block size")
    parser.add_argument('-tours', type=int, default=20, help="Max BKZ tours")
    parser.add_argument('-ghfactor', type=float, default=1, help="BKZ GH factor")
    parser.add_argument('-subtree-root-level', type=int, help="subtrees are rooted at this level")
    parser.add_argument('-tour', default=None)
    parser.add_argument('-block', default=None)
    args = parser.parse_args()

    data = json.load(open(args.subtree_stats_filename))
    treedata = None if not args.tree_stats_filename else json.load(open(args.tree_stats_filename))
    lwe = { 'n': args.n, 'q': args.q, 'sd': args.sd, 'm': args.m }
    bkz = { 'beta': args.beta, 'tours': args.tours, 'ghfactor': args.ghfactor }

    plots_path = f"{args.plot_dir}/{args.subtree_root_level}"
    os.system('mkdir -p '+ plots_path)
    gen_plots(data, plots_path, lwe, bkz, args.subtree_root_level, _tour=args.tour, _block=args.block, treedata=treedata)

def sims(lwe, bkz):
    n, q, sd, m = lwe['n'], lwe['q'], lwe['sd'], lwe['m']
    beta, tours, ghfactor = bkz['beta'], bkz['tours'], bkz['ghfactor']
    print(n, q, sd, m, beta, tours, ghfactor)

    lll_prof = LLLProfile(n, q, m)
    sim = Sim(lll_prof)
    save(sim.plot(), "test_lll.png")
    sim(beta, tours)
    save(sim.plot(), "test_bkz.png")
