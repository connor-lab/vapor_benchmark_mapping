import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
sns.set_style("ticks")
import sys
import pandas as pd
import matplotlib.gridspec as gridspec

def plot_hists(dfs):
    fig, axes = plt.subplots(nrows=4, ncols=3, sharey=True, sharex=True, figsize=(10,10))
    i = 0
    flat = axes.flatten()
    dfc = 0
    bins = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
    # Plot each of the columns (4,5,6,7) over nreads (3)
    for j in range(5,9):
        for dfc in range(len(dfs)):
            df = dfs[dfc]
            df.fillna
            arr = df.values
            axi = flat[i]
            # Exclude any datapoints that have zero reads, throw a warning that they are being excluded or we will divide by zero
            # If there are zero reads, the pipeline broke for that datapoint
            valind = [vi for vi in range(len(arr)) if arr[vi][3] != 0]
            print(j, dfc)
            print(len(arr)-len(valind), "datapoint(s) discarded due to zero reads")
            n_reads = arr[:,4][valind]
            n_mapped = arr[:,j][valind]
            perc_mapped = (n_mapped/n_reads).astype(float)
            print("Mean:", np.nanmean(perc_mapped))
            print("Std:", np.nanstd(perc_mapped))
#            axi.hist(perc_mapped,bins=10, density=True, histtype='step', linewidth=2.)
            axi.hist(perc_mapped,bins=bins, density=True)
            i += 1

    # Formatting
    cols = ['{}'.format(col) for col in ["Human", "Swine", "Avian"]]
    rows = ['{}'.format(row) for row in ['MINIMAP2', 'NGM', 'BWA', 'HISAT2']]
    for ax, col in zip(axes[0], cols):  
        ax.set_title(col)

    for ax, row in zip(axes[:,-1], rows):
        ax2 = ax.twinx()
        ax2.grid(False)
        ax2.get_yaxis().set_ticks([])
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax2.set_ylabel(row, labelpad=7)

    lastax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07)
    lastax.set_ylabel("Density",labelpad=3, size=14)
    lastax.set_xlabel("Proportion of Reads Mapping",labelpad=3, size=14)
    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    plt.savefig("mapping_tool_percs_h1n1_cali_huavsw.pdf", format="pdf")
    plt.show()

dfs = []
for fname in sys.argv[1:]:
    df = pd.read_csv(fname)
    dfs.append(df)

plot_hists(dfs)
