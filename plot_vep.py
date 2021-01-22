#!/usr/bin/python

#%%
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scipy.stats as stats
import statsmodels.stats.multicomp as mc


#%%
vep_table = pd.read_table('/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/vep_summary_stats.txt', sep = ',')
# %%
vep_table.head()

# %%
vep_table['miss_syn'] = vep_table['missense']/vep_table['syn']
vep_table['del_syn'] = vep_table['del']/vep_table['syn']
vep_table['LOF_syn'] = vep_table['stop_gained']/vep_table['syn']

#%%
vep_table.head()

#%%
min(vep_table['del_syn'])

#%%
def plot_mutationload(mutationload, label):
    order = ["F.Pyrenees", "F.Alps", "Iceland", "Sweden", "E.Greenland", "W.Greenland", "Svalbard", "Newfoundland", "England", "Norway"]
    ax = sns.boxplot(y = vep_table[mutationload], x = vep_table['pop'], width = 1, palette="colorblind", order = order)
    ax = sns.swarmplot(y = vep_table[mutationload], x = vep_table['pop'], size = 3, edgecolor='black', linewidth=0.9, order = order)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45)
    ax.set_ylabel(label)
    ax.set_xlabel('')
    ax.figure.savefig(f'{mutationload}.pdf', bbox_inches = 'tight')

#%%
plot_mutationload('miss_syn', label = 'Missense/Synonymous')
#%%
plot_mutationload('del_syn', label = 'Deleterious/Synonymous')
#%%
plot_mutationload('LOF_syn', label = 'Loss of function/Synonymous')

#%%
plot_mutationload('syn', label = 'Synonymous')
#%%
plot_mutationload('del', label = 'Deleterious')
#%%
plot_mutationload('stop_gained', label = 'Loss of function')


# %% #### Missense
model_miss = ols('miss_syn ~ C(Pop)', data = vep_table).fit()
aov_table_miss = sm.stats.anova_lm(model_miss, typ=2)
aov_table_miss

# %% #### Deleterious
model_del = ols('del_syn ~ C(Pop)', data = vep_table).fit()
aov_table_del = sm.stats.anova_lm(model_del, typ=2)
aov_table_del

# %% #### LOF
model_lof = ols('LOF_syn ~ C(Pop)', data = vep_table).fit()
aov_table_lof = sm.stats.anova_lm(model_lof, typ=2)
aov_table_lof

# %%
#Calculate effect sizes
# etasqe is the anova equivalent of R2. They are both effect sizes
"""
The function below was created specifically for the one-way ANOVA table results returned for Type II sum of squares
"""

def anova_table(aov):
    aov['mean_sq'] = aov[:]['sum_sq']/aov[:]['df']

    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])

    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*aov['mean_sq'][-1]))/(sum(aov['sum_sq'])+aov['mean_sq'][-1])

    cols = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[cols]
    return aov

#%%
anova_table(aov_table_miss)
anova_table(aov_table_del)
anova_table(aov_table_lof)

#### Checking assumptions
##Independence - Yes

#%%
##Normality of residuals using a shapiro-wilks test, however very conservative. NS = normal distributed residuals
stats.shapiro(model_miss.resid)
#%%
stats.shapiro(model_del.resid)
#%%
#%%
stats.shapiro(model_lof.resid)

# %%
### Another less conservative approach is visual inspection
def plot_normality(model, mutationload):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    normality_plot, stat = stats.probplot(model.resid, plot = plt, rvalue=True)
    ax.set_title("Prob plot of model resids "+ mutationload, fontsize = 20)
    ax.set

#%%
plot_normality(model_miss, mutationload = 'miss/syn')
#%%
plot_normality(model_del, mutationload = 'del/syn')
#%%
plot_normality(model_lof, mutationload = 'lof/syn')

# %%
## Homogeneity of variance using a levenes test, ns = homogenous variances

def stats_levene(mutationload):
    return stats.levene(vep_table[mutationload][vep_table['Pop'] == 'F.Pyrenees'],
    vep_table[mutationload][vep_table['Pop'] == 'F.Alps'],
    vep_table[mutationload][vep_table['Pop'] == 'Iceland'],
    vep_table[mutationload][vep_table['Pop'] == 'Sweden'],
    vep_table[mutationload][vep_table['Pop'] == 'E.Greenland'],
    vep_table[mutationload][vep_table['Pop'] == 'W.Greenland'],
    vep_table[mutationload][vep_table['Pop'] == 'Svalbard'],
    vep_table[mutationload][vep_table['Pop'] == 'Newfoundland'],
    vep_table[mutationload][vep_table['Pop'] == 'England'],
    vep_table[mutationload][vep_table['Pop'] == 'Norway'])

#%%
stats_levene('miss_syn')
#%%
stats_levene('del_syn')
#%%
stats_levene('LOF_syn')
#%%
#### Post-hoc testing to compare individual group means using Tukey honestly sign. difference (HSD)

def post_hoc_tukey(mutationload, label):
    comp = mc.MultiComparison(vep_table[mutationload], vep_table['Pop'])
    post_hoc_res = comp.tukeyhsd()
    plot = post_hoc_res.plot_simultaneous(ylabel=label, xlabel="Population")
    return post_hoc_res.summary(), plot
# %%
post_hoc_res_miss, plot_miss = post_hoc_tukey('miss_syn', label = 'miss/syn')
# %%
post_hoc_res_miss
# %%
post_hoc_res_del, plot_miss = post_hoc_tukey('del_syn', label = 'del/syn')
# %%
post_hoc_res_del
# %%
post_hoc_res_lof, plot_lof = post_hoc_tukey('LOF_syn', label = 'lof/syn')
# %%
post_hoc_res_lof
