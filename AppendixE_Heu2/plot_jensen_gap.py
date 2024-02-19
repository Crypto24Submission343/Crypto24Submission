import gc
import os 
import json
import profile
from random import gauss
from matplotlib import legend_handler
from sage.all import line, save, sqrt, log, prod, beta, list_plot, multi_graphics, histogram, scatter_plot, randint, seed
from numpy import mean
from bkz_simulators.sim import Sim, LLLProfile
import sys
sys.path.append("../LatticeHeuristics")
import enumTrees

def most_common(lst):
    cmn = max(set(lst), key=lst.count)
    qty = sum(int(_ == cmn) for _ in lst)
    return cmn, qty

def gaussian_heuristic(v, n, verbose=True):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """
    from math import pi, sqrt, e
    return ((pi * n)**(1./(2.*n))) * sqrt(n/(2 * pi * e)) * (v ** (1./n))
    from sage.all import RR
    v = RR(v)
    n = RR(n)
    a = ((pi * n)**(1./(2.*n))) * sqrt(n/(2 * pi * e))
    b = (v ** (1./n))
    if verbose:
        print("a", a)
        print("b", b)
        print("ab", a*b)
        print("og", ((pi * n)**(1./(2.*n))) * sqrt(n/(2 * pi * e)) * (v ** (1./n)))
    return a*b


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


def gen_plots(data, plots_dir, lwe, bkz, subtree_root_level, _tour=None, _block=None, treedata=None, dpi=150):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    sim = Sim(LLLProfile(n, q, m))

    plot_histograms = True

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
        subtree_avg = []
        sqrt_subtree_avg = []
        biased_sample_variance_subtree = []
        unbiased_sample_variance_subtree = []
        subtree_sizes_as_visited = []
        H_n_minus_k = {}
        theory_subtree_avg_sim = {}
        theory_subtree_avg_gs2 = {}
        for index in tour_stats_tuple:
            # if int(index) != 100:
            #     continue

            print(f"tour: {tour}, index: {index}, beta: {bkz['beta']}")

            sim_block_prof = sim.profile[int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
            sim_gh = profile_gh(sim_block_prof)
            if treedata:
                gs2 = treedata[tour][index]["gs"][int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
            sim_R = sqrt(bkz['ghfactor']) * sim_gh
            stats = tour_stats[index]
            avg = stats["avg"][subtree_root_level:]

            subtree_avg.append((int(index.split("-")[0]), stats["avg_tot"]))
            sqrt_subtree_avg.append((int(index.split("-")[0]), stats["avg_sqrt_tot"]))
            _b_samp_var = stats["avg_sqr_tot"] - stats["avg_tot"]**2
            biased_sample_variance_subtree.append((int(index.split("-")[0]), _b_samp_var))
            unbiased_sample_variance_subtree.append((int(index.split("-")[0]), stats["unb_samp_var_tot"]))
            try:
                subtree_sizes_as_visited.append((index, stats["subtree_sizes_as_visited"]))
                H_n_minus_k[index] = stats["H_n_minus_k"]
            except:
                plot_histograms = False

            from sage.all import pi, gamma
            k = subtree_root_level
            n = bkz['beta']
            def Ckh(k, h, sim):
                # formula obtained by direct integration.
                res = pi**(h/2) * sim_R**h * gamma(k/2+1)/gamma((k+h)/2+1) ########### dimension k? not n-k???????
                if sim:
                    covol = prod([sim_block_prof[i-1] for i in range(bkz['beta']-k, bkz['beta']-k-h, -1)])
                else:
                    covol = prod([gs2[i-1] for i in range(bkz['beta']-k, bkz['beta']-k-h, -1)])
                covol = sqrt(covol)
                res = res / covol
                return res

            def Nkh(k, sim):
                return float(sum(Ckh(k, i+1, sim) for i in range(len(avg))))

            try:
                theory_subtree_avg_sim[int(index.split('-')[0])] = Nkh(k, True)
                theory_subtree_avg_gs2[int(index.split('-')[0])] = Nkh(k, False)
            except:
                pass
            # print(theory_subtree_avg_sim)

        # """

            # print(index, k)
        g = line([])

        g += list_plot([(i, sqrt(j)) for i, j in subtree_avg],  marker="D", size=20, color="blue", legend_label=r"measured $\sqrt{\mathsf{avg}({\mathcal{T}(g \in Z_k)})}$", axes_labels=["block index $i$", ""], xmin=0, xmax=len(sim.profile)-k)
        g += list_plot([(i, j) for i, j in sqrt_subtree_avg], marker="D", size=20, color="green", legend_label=r"measured $\mathsf{avg}({\sqrt{\mathcal{T}(g \in Z_k)}})$", axes_labels=[ "block index $i$", ""], xmin=0, xmax=len(sim.profile)-k)
        # g += list_plot([(i, j) for i, j in biased_sample_variance_subtree], marker="D", size=20, color="red", legend_label=r"measured $\sqrt[4]{\frac{n-1}{n}  s^2({\mathcal{T}(g \in Z_k)})}$", axes_labels=[ "block index $i$", ""], xmin=0, xmax=len(sim.profile)-k)

        g.set_legend_options(loc='lower left')

        filename = f"{plots_dir}/{tour}"
        print(filename)
        # save(g, filename+".png", dpi=dpi)
        save(g, filename+".pdf")

        g += list_plot([(i, sqrt(sqrt(j))) for i, j in unbiased_sample_variance_subtree], marker="D", size=20, color="black", legend_label=r"measured $\sqrt[4]{s^2({\mathcal{T}(g \in Z_k)})}$", axes_labels=[ "block index $i$", ""], xmin=0, xmax=len(sim.profile)-k)
        filename = f"{plots_dir}/{tour}-incl-sigma"
        print(filename)
        # save(g, filename+".png", dpi=dpi)
        save(g, filename+".pdf")

        # """


        if plot_histograms and tour in [4, 9, 14]:
            for ind, subtree_sizes in subtree_sizes_as_visited:
                if int(ind.split("-")[0]) % 10 != 0 or len(subtree_sizes) <= 1:
                    continue
                _subtree_sizes = [_ + 1 for _ in subtree_sizes] # +1 due to root node
                gc.collect()

                # computing expected tree size
                _k = subtree_root_level
                _n = bkz['beta']
                # NkD = 2**float(enumTrees.NoPruning.log_avg_subtree_N_kD(_k, _n-_k, _n))
                # k+1 + ... + n
                # n - (k+1)            n - n 
                _k = subtree_root_level
                _n = bkz['beta']
                H_k = H_n_minus_k[ind][::-1]
                print(H_k)
                plot_expectation = True
                try:
                    NkD = sum(H_k[i]/H_k[_k-1] for i in range(_k, _n))
                except:
                    plot_expectation = False
                if plot_expectation:
                    print("################", int(ind.split('-')[0]), NkD)

                # sizes as visited
                filename = f"{plots_dir}/visit/{tour}-{ind}-subtree-sizes-as-visited"
                g = line([])
                # g += line([(0, _lb_sT_size_sim), (len(_subtree_sizes), _lb_sT_size_sim)], color="red", legend_label="sim lb", thickness=2)
                # g += line([(0, _lb_sT_size_gs2), (len(_subtree_sizes), _lb_sT_size_gs2)], color="black", legend_label="gs2 lb", thickness=2)
                g += line([(0, mean(_subtree_sizes)), (len(_subtree_sizes), mean(_subtree_sizes))], color="green", legend_label="average", thickness=1, linestyle="solid")
                if plot_expectation:
                    g += line([(0, NkD), (len(_subtree_sizes), NkD)], color="black", legend_label="$\sum{H_{k+h}/{H_k}}$", thickness=1, linestyle="--")
                g += scatter_plot([(i, _subtree_sizes[i]) for i in range(len(_subtree_sizes))], dpi=dpi, marker='.', alpha=.5)
                # for i in range(1, 20):
                #     g += line([(0, i), (500, i)])
                print(filename)
                # save(g, filename+".png", dpi=dpi)
                save(g, filename+".pdf")


                # scatter of random subsets and of serial subsets
                                # sizes as visited
                filename = f"{plots_dir}/comparison/{tour}-{ind}-virtual-subtrees"
                # collect sequential sets of 2^y nodes, except for first subtree
                virtual_trees_size = {4: [], 6: [], 8: []}
                for y in virtual_trees_size:
                    for i in range(1, len(_subtree_sizes), 2**y):
                        virtual_tree = _subtree_sizes[i:] if i+2**y >= len(_subtree_sizes) else _subtree_sizes[i:i+2**y]
                        virtual_trees_size[y].append(sum(virtual_tree))

                    g = line([])
                    g += scatter_plot([(x, virtual_trees_size[y][x]) for x in range(len(virtual_trees_size[y]))], marker='.')
                    g += line([(0, mean(virtual_trees_size[y])), (len(virtual_trees_size[y]), mean(virtual_trees_size[y]))], color="green", legend_label="mean", thickness=1)
                    if plot_expectation:
                        g += line([(0, 2**y * NkD), (len(virtual_trees_size[y]), 2**y * NkD)], color="black", legend_label="$2^y \cdot \sum{H_{k+h}/{H_k}}$", thickness=1, linestyle="--")
                    save(g, f"{filename}-y{y}.pdf")


                # full histogram, except for first, largest subtree
                filename = f"{plots_dir}/histogram/{tour}-{ind}-histogram-subtree-sizes"
                g = line([])
                g += histogram(_subtree_sizes, bins=int(max(_subtree_sizes)+1), dpi=dpi)
                print(_subtree_sizes)
                print(most_common(_subtree_sizes))
                print(H_n_minus_k[index])
                print()
                _mean = mean(_subtree_sizes)
                _height = most_common(_subtree_sizes)[1]
                g += line([(mean(_subtree_sizes), 0), (mean(_subtree_sizes), _height)], color="green", legend_label="mean",)
                if plot_expectation:
                    g += line([(NkD, 0), (NkD, _height)], color="black", legend_label="$\sum{H_{k+h}/{H_k}}$", linestyle="--")

                print(filename)
                # save(g, filename+".png", dpi=dpi)
                save(g, filename+".pdf")

                # histogram except first very large subtree
                filename = f"{plots_dir}/histogram-except-first/{tour}-{ind}-histogram-subtree-sizes"
                g = histogram(_subtree_sizes[1:], bins=int(max(_subtree_sizes[1:])+1), dpi=dpi)
                g += line([(mean(_subtree_sizes), 0), (mean(_subtree_sizes), _height)], color="green", legend_label="mean",)
                if plot_expectation:
                    g += line([(NkD, 0), (NkD, _height)], color="black", legend_label="$\sum{H_{k+h}/{H_k}}$", linestyle="--")
                print(filename)
                # save(g, filename+".png", dpi=dpi)
                save(g, filename+".pdf")

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
    os.system('mkdir -p '+ plots_path + '/histogram')
    os.system('mkdir -p '+ plots_path + '/histogram-except-first')
    os.system('mkdir -p '+ plots_path + '/visit')
    os.system('mkdir -p '+ plots_path + '/comparison')
    with seed(0xdeadbeef):
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
