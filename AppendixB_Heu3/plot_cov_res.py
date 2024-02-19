import numpy
import math
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
import os 
import statistics as st

ss_dim = 1
ss_k = 1 
ss_kh = 1

def read_input(fn):
    """
    Read cov_res.txt file and load data into dictionary.

    [dim][k][k+h] --> [
                        Cov(Zk^-1,Zkh) ,
                        Cov(Zk^-1,Zkh)*(H_{k+h}*H_k^-1)^-1 ,
                        (E[X]E[Y]^-1)E[XdivY]^-1 ,
                    ]    
    :param      fn:   The function
    :type       fn:   Function
    """

    file = open(fn, 'r')

    lines = (line.rstrip() for line in file)  
    lines = list(line for line in lines if line) # Non-blank lines in a list
    
    dim_n = 0
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(None)))

    def parse_line(line, dim_n):

        cols = line.split(';')
        k = int(cols[0].split(':')[1].strip())
        k_plus_h = int(cols[1].split(':')[1].strip())

        cov_Zk_Zkh = float(cols[2].split(':')[1].strip())
        cov_Zk_Zkh_Hkh_Hk = float(cols[3].split(':')[1].strip()) 
        Ex_Ey = float(cols[4].split(':')[1].strip()) 
        #
        H_kh = float(cols[5].split(':')[1].strip())
        H_k = float(cols[6].split(':')[1].strip()) 
        Exy = float(cols[7].split(':')[1].strip())  

        data[dim_n][k][k_plus_h] = [cov_Zk_Zkh, cov_Zk_Zkh_Hkh_Hk, Ex_Ey, H_kh, H_k, Exy]


    data_params = {}
    for line in lines:
        if line.startswith('tries'):
            parameter = line.split(',')
            p_string = ''
            for p in parameter:
                p_string += p.split('=')[0].strip()
                p_string += p.split('=')[1] + '-'
            dim_n = int(parameter[1].split('=')[1])

            data_params[dim_n] = p_string[:-1]

        if line.startswith('k'):
            parse_line(line, dim_n)
    return data, data_params


def setup_plot(dim, fn_path, fn_out='', yscale='symlog', ylabel='Value'):
    full_path = fn_path
    os.system('mkdir -p ' + f'{full_path}')
    os.system('touch ' + f'{full_path}' + '/PLACEHOLDER')

    fig, ax = plt.subplots()
    ax.set_yscale(yscale)

    return fig, ax, full_path

def write_plot(ax, fig, label, full_path, fn_out, dim, clear=True, add_span = True):

    if label == 'cov':
        min_y = max(ax.get_ylim()[0]-0.05, -1)
        ax.axhspan(ymin=min_y, ymax=0, alpha=.05, label="[-1, 0]", color='red')
        plt.ylabel('$|Cov(1/|Z_k |, |Z_{k+h}|)|$ over $h \in [k+1, n]$')


    if label == 'exp':
        #min_y = min(ax.get_ylim()[0]-0.05, -.5)
        #ax.axhspan(ymin=min_y, ymax=1.5, alpha=.05, label="[$\infty$, 1]", color='red')
        #plt.ylim(min_y, ax.get_ylim()[1])
        plt.ylim(.4, 1.6)

        ax.set_yscale('linear')
        plt.ylabel('$\\frac{H_{k+h}}{H_{k}}\ /\ S_{k,h}$')
    
    #ax.legend(loc='upper right')#, bbox_to_anchor=(1.04, 0.5))

    plt.savefig(f'{full_path}/' + str(fn_out), dpi=300, transparent=False, bbox_inches='tight')
    print(f"Wrote file " + f'{full_path}/' + str(fn_out))

    # if label == 'cov':
    #     plt.ylim(-2, .5)
    #     ax.set_yticks(np.arange(-2, .5, step=0.5))
    #     ax.set_yticklabels(np.arange(-2, .5, step=.5))
    #     #
    #     plt.xlim(math.floor(dim/2)-1, dim)
    #     #
    #     plt.savefig(f'{full_path}/' + 'zoom_' + str(fn_out), dpi=300, transparent=False, bbox_inches='tight')
    #     print(f"Wrote file " + f'{full_path}/' + str(fn_out))

    # if label == 'exp':
    #     plt.ylim(-.5, 1.5)
    #     ax.set_yticks(np.arange(-.5, 1.5, step=0.5))
    #     ax.set_yticklabels(np.arange(-.5, 1.5, step=.5))
    #     #
    #     plt.xlim(math.floor(dim/2)-1, dim)
    #     plt.savefig(f'{full_path}/' + 'zoom_' + str(fn_out), dpi=300, transparent=False, bbox_inches='tight')
    #     print(f"Wrote file " + f'{full_path}/' + str(fn_out))

    if clear:
        plt.close()
        plt.cla()
        plt.clf()
    


def iterate_data(data, col, ss_dim=4, ss_k=4, ss_kh=4):
    reduced_dim = list(data.keys())[::ss_dim]
    for dim in reduced_dim:

        # Root level
        reduced_k = list(data[dim].keys())[::ss_k]
        for k in reduced_k:

            # Height level
            reduced_k_plus_h = list(data[dim][k].keys())[::ss_kh]
            for k_plus_h in reduced_k_plus_h:
                yield dim, k, k_plus_h, data[dim][k][k_plus_h][col]


def plot_group_kh(dim_k_kh_values, col, fn_out, data_params, label):
    """
    Plot by grouping all values k+1, k+2, ..., k+h, ..., n into one boxplot
    
    :param      dim_k_kh_values:  The dim k kh values
    :type       dim_k_kh_values:  { type_description }
    :param      col:              The col
    :type       col:              { type_description }
    :param      fn_out:           The function out
    :type       fn_out:           { type_description }
    """

    dim_to_level_to_value = defaultdict(lambda: defaultdict(lambda: []))

    reduced_dim = list(data.keys())[::ss_dim]
    for dim in reduced_dim:

        # Root level
        reduced_k = list(data[dim].keys())[::ss_k]
        for k in reduced_k:

            # Height level
            reduced_k_plus_h = list(data[dim][k].keys())[::ss_kh]
            for k_plus_h in reduced_k_plus_h:
                if label == 'exp':
                    EXEY_EXY1 = (data[dim][k][k_plus_h][3] / data[dim][k][k_plus_h][4]) / (data[dim][k][k_plus_h][5])  # + 1
                    dim_to_level_to_value[dim][k].append(EXEY_EXY1)
                else:
                    dim_to_level_to_value[dim][k].append(data[dim][k][k_plus_h][col])

    # Set x-positions for boxes
    x_pos_range = np.arange(len(dim_to_level_to_value)) / (max(len(dim_to_level_to_value),2) - 1)
    x_pos = (x_pos_range * .75) + 0.75

    labeled = False
    for i, (dim, k_to_v) in enumerate(dim_to_level_to_value.items()):
        # No pruning
        fig, ax, full_path = setup_plot(dim, 'group_forall_kh')

        # Linear pruning
        # fig, ax, full_path = setup_plot(dim, 'group_forall_kh_lin')

        palette = sns.color_palette(None, max(dim_to_level_to_value[dim].keys()))   

        ax.grid(which='major', alpha=0.5)
        plt.xlabel("k")
        #plt.title("For each $k$: $Cov(1/|Z_k|, |Z_{k+h}|)$ over $h \in [k+1, n]$")

        # Plot
        pos = [x_pos[i] + j * ss_k for j in range(len(k_to_v.keys()))]

        bp = plt.boxplot(k_to_v.values(), #widths=0.6 / len(k_to_v), 
            patch_artist=True, # Required for coloring 
            boxprops=dict(facecolor=palette[i], color=palette[i]), # Sets color 
            meanline=True, showmeans=True, 
            flierprops={'markersize': 1}, 
            #labels=list(k_to_v.keys()),
            #labels=[' '] * len(h_to_v.keys()),
            positions=pos
        )

        ax.set_xticks(np.arange(0, dim, step=10))
        ax.set_xticklabels(np.arange(0, dim, step=10))
        #ax.xaxis.set_major_locator(MaxNLocator(integer=True))


        # bp = plt.violinplot(k_to_v.values(), #widths=0.6 / len(k_to_v), 
        #     #patch_artist=True, # Required for coloring 
        #     #boxprops=dict(facecolor=palette[i], color=palette[i]), # Sets color 
        #     #
        #     # manage_ticks=False, meanline=True, 
        #     #labels=list(k_to_v.keys()),
        #     #labels=[' '] * len(h_to_v.keys()),
        #     positions=pos
        # )
        
        # Mark k >= n/2
        ax.axvline(x=dim/2, dashes=[4, 8], color=palette[i], alpha=1, lw=1, label='k $\geq$ n/2')

        write_plot(ax, fig, label, full_path, data_params[dim] + '_' + fn_out, dim=dim)

def plot_scatter_k_kh(dim_k_kh_values, col, fn_out, data_params, label):

    dim_to_level_to_value = defaultdict(lambda: defaultdict())

    reduced_dim = list(data.keys())[::ss_dim]
    for dim in reduced_dim:
        # Root level
        reduced_k = list(data[dim].keys())[::ss_k]
        for k in reduced_k:

            xs = [] 
            ys = []

            # Height level
            reduced_k_plus_h = list(data[dim][k].keys())[::ss_kh]
            for k_plus_h in reduced_k_plus_h:
                xs.append(str(k) + '-' + str(k_plus_h))

                if label == 'exp':
                    EXEY_EXY1 = (data[dim][k][k_plus_h][3] / data[dim][k][k_plus_h][4]) / (data[dim][k][k_plus_h][5]) # + 1
                    ys.append(EXEY_EXY1)
                else:
                    ys.append(data[dim][k][k_plus_h][col])

            dim_to_level_to_value[dim][k] = [xs, ys]
    # End Bucket Loop

    # Generate Colors
    palette = sns.color_palette(None, 50)    


    labeled = False
    for i, (dim, k_to_xy) in enumerate(dim_to_level_to_value.items()):

        fig, ax, full_path,  = setup_plot(dim, 'scatter_kh')

        plt.xlabel("k,k+h")
        #plt.title(f"One point for every tuple (k, k+h), one color per k")

        for k, xy in k_to_xy.items():
            ax.scatter(xy[0], xy[1], s=0.75)

        write_plot(ax, fig, label, full_path, data_params[dim] + '_' + fn_out, dim=dim)

def plot_lines_from_k(dim_k_kh_values, col, fn_out, data_params, label, dir_out):
    dim_to_level_to_value = defaultdict(lambda: defaultdict(lambda: []))

    reduced_dim = list(data.keys())[::ss_dim]
    for dim in reduced_dim:
        
        # Root level
        reduced_k = list(data[dim].keys())[::ss_k]
        for k in reduced_k:
            xs = [] 
            ys = []

            # Height level
            reduced_k_plus_h = list(data[dim][k].keys())[::ss_kh]
            for k_plus_h in reduced_k_plus_h:
                xs.append(k_plus_h)

                if label == 'exp':
                    EXEY_EXY1 = data[dim][k][k_plus_h][2]
                    #EXEY_EXY1 = (data[dim][k][k_plus_h][3] / data[dim][k][k_plus_h][4]) / (data[dim][k][k_plus_h][5] + 1)
                    ys.append(EXEY_EXY1)
                else:
                    ys.append(data[dim][k][k_plus_h][col])
            dim_to_level_to_value[dim][k] = [xs, ys]
    # End Bucket Loop

    # Generate Colors
    palette = sns.color_palette(None, max(dim_k_kh_values.keys()))    


    labeled = False
    for i, (dim, k_to_xy) in enumerate(dim_to_level_to_value.items()):

        # No Pruning
        fig, ax, full_path = setup_plot(dim, f'{dir_out}/plot_lines_kh', fn_out)

        # Linear Pruning
        # fig, ax, full_path = setup_plot(dim, 'plot_lines_kh_lin', fn_out)

        ax.grid(which='major', alpha=0.5)
        plt.xlabel("k")
        #plt.title("For each $k$: $Cov(1/|Z_k|, |Z_{k+h}|)$ over $h \in [k+1, n]$")

        for i, (k,xy) in enumerate(k_to_xy.items()): 
            alpha = 1
            color = palette[i]
            # if k <= dim/2:
            #     alpha = 0.4
            #     if label == 'exp':
            #         alpha = 0.2
            #     color = 'grey'

            ax.plot(xy[0], xy[1], color=color, linewidth=.5, alpha=alpha)

            plt.plot(k+1, xy[1][0], marker='o', color=color, markersize=2.) 
            #plt.text(x=xy[0][0]-1, y=xy[1][0]-.04, s=f"k={k}", fontsize=4) #xy[0]

            # anot = ax.annotate(f"k={k}", xy=(xy[0][0], xy[1][0]), xytext=(xy[0][0]-1, xy[1][0]-.1),
            #                color='black', ha='center', va='center', fontsize=4, alpha=1,
            #                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", linewidth=0.5))

            # xy=(x_coord, y_coord), xytext=(x_coord + x_offset, y_coord + y_offset),
        # # Mark k >= n/2
        #ax.axvline(x=dim/2, dashes=[4, 8], color=palette[i], alpha=1, lw=1, label='k $\geq$ n/2')

        write_plot(ax, fig, label, full_path, data_params[dim] + '_' + fn_out, dim=dim)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default="cov_res.txt")
    parser.add_argument('-o', '--outputdir', default=".")
    parser.add_argument('-e', '--ext', default="pdf")
    args = parser.parse_args()

    data, data_params = read_input(args.input)
    #plot_group_kh(data, col=0, fn_out='covariance.pdf', data_params=data_params, label='cov')
    #plot_group_kh(data, col=2, fn_out='expectedVal.pdf', data_params=data_params, label='exp')
    
    #plot_scatter_k_kh(data, col=0, fn_out='Cov(Zk^-1,Zkh).pdf', data_params=data_params, label='cov')
    #plot_scatter_k_kh(data, col=2, fn_out='(E[X]E[Y]^-1)E[XdivY]^-1.pdf', data_params=data_params, label='exp')

    #plot_lines_from_k(data, col=0, fn_out='covariance.pdf', data_params=data_params, label='cov') #Cov(Zk^-1,Zkh).pdf
    plot_lines_from_k(data, col=2, fn_out=f'expectedVal.{args.ext}', data_params=data_params, label='exp', dir_out=args.outputdir) #(E[X]E[Y]^-1)E[XdivY]^-1.pdf

    
