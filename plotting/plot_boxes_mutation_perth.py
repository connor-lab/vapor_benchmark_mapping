import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
import sys
import pandas as pd
import matplotlib.gridspec as gridspec
sns.set_style("ticks")

def plot_boxes(df, titles="MINIMAP2,NGM,BWA,HISAT2".split(",")):
    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True, figsize=(10,10))
    i = 0
    flat = axes.flatten()
    dfc = 0
    df.fillna(0, inplace=True)
    for j in range(4,8):
        arr = df.values
        axi = flat[i]
        nmapped = arr[:,j]
        ncol = arr[:,3]
        mutcol = arr[:,0]/100.

        sns.boxplot(x=mutcol, y=100.*nmapped/ncol, ax=axi)

        
        axi.set_title(titles[i])
        i += 1

    lastax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07)
    lastax.set_ylabel("Percentage Reads Recovered",labelpad=10,size=14)
    lastax.set_xlabel("Mutation Probability",labelpad=10,size=14)
    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    plt.savefig("boxes_perth_mutation.pdf", format="pdf", dpi=900)
#    plt.show()

df = pd.read_csv(sys.argv[1])
plot_boxes(df)
