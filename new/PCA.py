#!/usr/bin/python

#%%
import random
random.seed(42)
import time
import numpy as np
np.random.seed(42)
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas
import allel; print('scikit-allel', allel.__version__)

# %%
muta = allel.read_vcf('LL_SNPs_final_autosomes.recode.vcf', fields = '*')

# %%
sorted(muta.keys())

# %%
muta['calldata/GT'].shape[0]

# %%
gt = muta['calldata/GT']
gt = allel.GenotypeArray(gt)
len(gt)

# %%
gn = gt.to_n_alt()
gn

# %%
coords1, model1 = allel.pca(gn, n_components=10, scaler=None)


# %%
df_samples = pandas.read_csv('LL_pop.txt', delimiter='\t', header = None)
df_samples.head()


# %%
populations = df_samples.iloc[:, 1].unique()
len(populations)
	
# %%
pop_colours = {
    'French.alps': '#FF0000',
    'E.Greenland': '#008000',
    'Iceland': '#00FFFF',
    'W.Greenland': '#90EE90',
    'Sweden(Jamtland)': '#FFA500',
    'Pyrenees': '#8B0000',
    'Svalbard': '#1E90FF'
}
# %%
pop_colours = {
    'Norway': '#FF0000',
    'Newfoundland': '#008000',
    'England(scoticus)': '#00FFFF'
}


#%%
def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))
    

#%%
def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = df_samples.iloc[:, 1].values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig




# %%
f = fig_pca(coords1, model1, 'PCA.')


# %%
f.savefig("LL_final.pdf")



# %%
