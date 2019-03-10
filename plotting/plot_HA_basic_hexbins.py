import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import matplotlib.gridspec as gridspec
import seaborn as sns; sns.set_style("ticks")

def plot_hexbins(df, fname, titles="MINIMAP2,NGM,BWA,HISAT2".split(",")):
    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True, figsize=(10,10))
    i = 0
    flat = axes.flatten()
    dfc = 0
    ims = []
    for j in range(6,10):
        arr = df.values
        axi = flat[i]
        col = arr[:,j]
        ncol = arr[:,5]
        idcol = arr[:,3]
        pidstr = idcol
        # Exclude any datapoints for which no reads were produced due to an error (otherwise we will divide by zero), report how many
        invalind = [i for i in range(len(pidstr)) if ncol[i] == 0]
        print("Excluded", len(invalind), "datapoints that had zero reads due to a broken pipeline")
        validind = [i for i in range(len(pidstr)) if pidstr[i] != "False" and ncol[i] != 0]
        pidstr = pidstr[validind]
        idcol = idcol[validind]
        ncol = ncol[validind]
        col = col[validind]
        reccol = 100.*col/ncol
        print("Mean:", titles[i], np.mean(reccol))
        im = axi.hexbin(idcol, reccol, cmap='inferno', gridsize=50, bins='log')
        axi.set_title(titles[i])
        i += 1

    # Format
    lastax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.subplots_adjust(left=0.07, right=0.9, top=0.95, bottom=0.07)
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.07, 0.03, 0.88])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel('Log Count', rotation=270, labelpad=13, size=14)
    lastax.set_ylabel("Percentage Reads Recovered",labelpad=12, size=14)
    lastax.set_xlabel("Percentage Identity",labelpad=10, size=14)
    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    plt.savefig(fname, format="pdf", dpi=900)

df = pd.read_csv(sys.argv[1])
fname = sys.argv[2]
plot_hexbins(df, fname)
