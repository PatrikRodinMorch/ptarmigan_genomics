#!/usr/bin/python

#%%
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Script to plot runs of homozygosity output from plink.


#%%
roh_table = pd.read_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/roh_table2.csv')
# %%
roh_table.head()
# %%

# getting genome size in to kb
autosomal_genome_size_kb = 960796788/1000
autosomal_genome_size_kb

#%%
roh_table['Froh_all'] = roh_table['KB_tot'] / autosomal_genome_size_kb
roh_table.head()
#%%
roh_table['Froh_1MB'] = roh_table['KB_long'] / autosomal_genome_size_kb
roh_table.head()
#%%
roh_table['Froh_100Kb'] = roh_table['KB_short'] / autosomal_genome_size_kb
roh_table.head()

# %%
def plot_roh(roh, label):
    order = ["F.Pyrenees", "F.Alps", "Iceland", "Sweden", "E.Greenland", "W.Greenland", "Svalbard", "Newfoundland", "England", "Norway"]
    ax = sns.boxplot(y = roh_table[roh], x = roh_table['Pop'], width = 1, palette="colorblind",  order = order)
    ax = sns.swarmplot(y = roh_table[roh], x = roh_table['Pop'], size = 3, edgecolor='black', linewidth=0.9, order = order)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45)
    ax.set_ylabel(label)
    ax.set_xlabel('')
    ax.figure.savefig(f'{roh}.pdf', bbox_inches = 'tight')

#%%
roh_table['KB_long_MB'] = roh_table['KB_long']/1000
plot_roh('KB_long_MB', label = 'Long ROH (Mb)')
#%%
roh_table['KB_short_MB'] = roh_table['KB_short']/1000
plot_roh('KB_short_MB', label = 'Short ROH (Mb)')

#%%
plot_roh('Froh_all', label = 'Froh')

#%%
plot_roh('Froh_1MB', label = 'Froh >1Mb')

#%%
plot_roh('Froh_100Kb', label = 'Froh >100Kb')

# %%
roh_table['KB_tot_MB'] = roh_table['KB_tot']/1000
plot_roh('KB_tot_MB', label = 'Total ROH length (Mb)')

# %%
def plot_roh_scatter(roh, label, roh_class):
    #order = ["F.Pyrenees", "F.Alps", "Iceland", "Sweden", "E.Greenland", "W.Greenland", "Svalbard", "Newfoundland", "England", "Norway"]
    ax = sns.scatterplot(y = roh_table[roh], x = roh_table[roh_class], palette="colorblind", hue = roh_table['Pop'], y_jitter= True, x_jitter= True)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1, fontsize = '10')
    ax.set_ylabel(label)
    ax.set_xlabel('Total number of ROH')
    ax.figure.savefig(f'{roh_class}.pdf', bbox_inches = 'tight')


#%%
plot_roh_scatter('KB_long_MB', 'ROH > 1Mb',  'ROH_number_long')

#%%
plot_roh_scatter('KB_short_MB', '100Kb < ROH < 1Mb',  'ROH_number_short')

#%%
plot_roh_scatter('KB_tot_MB', 'ROH > 100Kb',  'ROH_number_tot')
# %%
